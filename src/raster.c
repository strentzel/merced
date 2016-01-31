#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef NDEBUG
#include <math.h>
#endif

#include "dscan.h"
#include "dscanimp.h"

void rasterize( convexpolygon_t* poly, int* bitmap, int w, int h )
{
    int lowestLeftIndex = 0;
    int highestLeftIndex = 0;
    int lowestRightIndex = 0;
    int highestRightIndex = 0;
    int i, y;
    int yStart, yEnd;
    int left, nextLeft, right, nextRight;
    int nextMinX, nextMaxX;
    for( i = 0; i < poly->mPointsCnt; ++i )
    {
        float x = poly->mPoints[i].x;
        float y = poly->mPoints[i].y;
        float llx = poly->mPoints[lowestLeftIndex].x;
        float lly = poly->mPoints[lowestLeftIndex].y;
        float hlx = poly->mPoints[highestLeftIndex].x;
        float hly = poly->mPoints[highestLeftIndex].y;
        float lrx = poly->mPoints[lowestRightIndex].x;
        float lry = poly->mPoints[lowestRightIndex].y;
        float hrx = poly->mPoints[highestRightIndex].x;
        float hry = poly->mPoints[highestRightIndex].y;
        if( y < lly || (y == lly && x < llx) ) lowestLeftIndex = i;
        if( y > hly || (y == hly && x < hlx) ) highestLeftIndex = i;
        if( y < lry || (y == lry && x > lrx) ) lowestRightIndex = i;
        if( y > hry || (y == hry && x > hrx) ) highestRightIndex = i;
    }

    yStart = (int) poly->mPoints[lowestLeftIndex].y;
    yEnd = (int) poly->mPoints[highestRightIndex].y;

    left = lowestLeftIndex;
    nextLeft = (left + poly->mPointsCnt - 1) % poly->mPointsCnt;
    right = lowestRightIndex;
    nextRight = (right + 1) % poly->mPointsCnt;


    nextMinX = poly->mPoints[left].x;
    nextMaxX = poly->mPoints[right].x;

    if( yEnd >= h ) yEnd = h-1;

    for( y = yStart; y <= yEnd; ++y )
    {
        int currTop = y + 1;
        int minX = nextMinX;
        int maxX = nextMaxX;

        while( poly->mPoints[nextLeft].y < currTop &&
               nextLeft != highestLeftIndex)
        {
            int nextX = (int) poly->mPoints[nextLeft].x;
            if( nextX < minX ) minX = nextX;
            left = nextLeft;
            nextLeft = (nextLeft + poly->mPointsCnt - 1) % poly->mPointsCnt;
        }

        if( y != yEnd )
        {
            float y0 = currTop;
            float sx = poly->mPoints[left].x;
            float sy = poly->mPoints[left].y;
            float tx = poly->mPoints[nextLeft].x;
            float ty = poly->mPoints[nextLeft].y;
            /* Note: computation should be valid even for tiny ty-sy since 
             * that implies tiny y0-sy since y0-sy < ty-sy */
            assert( fabs(y0-sy) <= fabs(ty-sy) );
            nextMinX = ((y0 - sy) / (ty - sy)) * (tx - sx) + sx;
        }
                                                                              
        if( nextMinX < minX ) minX = nextMinX;

        while( poly->mPoints[nextRight].y < currTop && 
               nextRight != highestRightIndex)
        {
            int nextX = (int) poly->mPoints[nextRight].x;
            if( nextX > maxX ) maxX = nextX;
            right = nextRight;
            nextRight = (nextRight + 1) % poly->mPointsCnt;
        }

        if( y != yEnd )
        {
            float y0 = currTop;
            float sx = poly->mPoints[right].x;
            float sy = poly->mPoints[right].y;
            float tx = poly->mPoints[nextRight].x;
            float ty = poly->mPoints[nextRight].y;
            assert( fabs(y0-sy) <= fabs(ty-sy) );
            nextMaxX = ((y0 - sy) / (ty - sy)) * (tx - sx) + sx;
        }
                                                                              
        if( nextMaxX > maxX ) maxX = nextMaxX;
        
        if( y >= 0 )
        {
            int x;
            if( minX < 0 ) minX = 0;
            if( maxX >= w ) maxX = w-1;
        
            for( x = minX; x <= maxX; ++x )
            {
                int* buf;
                assert( x >= 0 && x < w && y >= 0 && y < h );
                buf = bitmap + y * w + x;
                *buf = 1;
            }
        }
    }
}

static void Fail()
{
    printf( "FAIL\n" );
    exit(-1);
}

typedef struct check_t
{
    int x;
    int y;
    int expected;
} check_t;

void test_rasterize( int count, int check, int dump, ... )
{
    va_list ap;
    vector2_t* points = malloc( sizeof(vector2_t) * count );
    check_t* checks = malloc( sizeof(check_t) * check );
    int i;
    int w = 0, h = 0;
    dscan_convexpolygon_t poly;
    int* data;
    va_start( ap, dump );
    for( i = 0; i < count; ++i )
    {
        points[i].x = va_arg( ap, double );
        points[i].y = va_arg( ap, double );
    }
    for( i = 0; i < check; ++i )
    {
        checks[i].x = va_arg( ap, int );
        checks[i].y = va_arg( ap, int );
        checks[i].expected = va_arg( ap, int );
    }
    va_end( ap );

    for( i = 0; i < count; ++i )
    {
        if( ceilf(points[i].x) > w ) w = (int) ceilf(points[i].x);
        if( ceilf(points[i].y) > h ) h = (int) ceilf(points[i].y);
    }
    w *= 2;
    h *= 2;

    poly = dscan_convexpolygon_create( points, count );
    data = malloc( sizeof(int) * w * h );
    for( i = 0; i < w*h; ++i ) data[i] = 0;

    rasterize( poly, data, w, h );
    for( i = 0; i < check; ++i )
    {
        if( checks[i].expected != data[checks[i].y*w + checks[i].x] ) 
        {
            Fail();
        }
    }

    if( dump )
    {
        int x, y;
        for( y = h-1; y >= 0; --y )
        {
            for( x = 0; x < w; ++x )
            {
                printf( "%d", data[y*w+x] );
            }
            printf( "\n" );
        }
    }
}

int main( int argc, char** argv )
{
    test_rasterize( 5, 0, 1,
                    2.0, 1.99999,
                    10.0, 1.99999,
                    15.0, 2.0,
                    10.0, 5.0,
                    2.0, 5.0 );

    test_rasterize( 4, 1, 0,
                    4.5, 6.5,
                    4.3, 6.4,
                    4.4, 6.18,
                    4.6, 6.6,
                    4, 6, 1 );

    test_rasterize( 4, 6, 0,
                    1.5, 18.5,
                    1.4, 10.0,
                    1.5, 1.5,
                    1.6, 10.0,
                    1, 1, 1,
                    1, 18, 1,
                    1, 0, 0,
                    1, 19, 0,
                    0, 10, 0,
                    2, 10, 0 );

    return 0;
}
