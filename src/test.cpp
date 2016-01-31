#include <fenv.h>
#include <float.h>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <list>
#include <set>
#include <vector>
using namespace std;

extern "C" 
{
#include "dscan.h"
#include "dscanimp.h"
}

namespace
{
    int XITERS = 1000;
    int CVXITERS = 1000;
}

bool gExitCalled;
extern "C" void catch_exit(int)
{
    gExitCalled = 1;
}

void Fail()
{
    cerr << "FAIL" << endl;
    exit(-1);
}

void Check( bool b )
{
    if( !b )
    {
        Fail();
    }
}

template<typename T>
void CheckEquiv( const T& s, const T& t )
{
    if( s < t || t < s )
    {
        Fail();
    }
}

void CheckEquiv( const vector3_t& s, const vector3_t& t )
{
    if( s.x != t.x || s.y != t.y || s.z != t.z )
    {
        Fail();
    }
}

void CheckEquiv( const vector2_t& s, const vector2_t& t )
{
    if( s.x != t.x || s.y != t.y )
    {
        Fail();
    }
}

template<typename T, typename DIFFT>
void CheckEquiv( const T& s, const T& t, DIFFT tolerance )
{
    if( std::abs( s - t ) > tolerance )
    {
        Fail();
    }
}

void CheckEquiv( const vector2_t& s, const vector2_t t, float tolerance )
{
    vector2_t dist = minus2( s, t );
    Check( norm2( dist ) < tolerance );
}

void CheckEquiv( const convexpolygon_t& s,
                 const convexpolygon_t& t, float tolerance )
{
    if( s.mPointsCnt != t.mPointsCnt ) Fail();

    for( int i = 0; i < s.mPointsCnt; ++i )
    {
        bool failed = 0;
        for( int j = 0, k = i; j < t.mPointsCnt; ++j, k=(k+1)%s.mPointsCnt )
        {
            vector2_t dist = minus2( s.mPoints[k], t.mPoints[j] );
            if( norm2( dist ) > tolerance )
            {
                failed = 1;
                break;
            }
        }
        if( !failed ) return;
    }
    if( s.mPointsCnt > 0 )
    {
        Fail();
    }
}

template<typename T>
void CheckEquiv( const set<T>& s, const set<T>& t )
{
    typedef typename set<T>::const_iterator iter_t;

    CheckEquiv( s.size(), t.size() );
    for( iter_t si = s.begin(), ti = t.begin(); si != s.end(); ++si, ++ti )
    {
        CheckEquiv( *si, *ti );
    }
}

vector3_t MakeVector( float x, float y, float z )
{
    vector3_t retval;
    retval.x = x;
    retval.y = y;
    retval.z = z;
    return retval;
}

matrix3_t MakeHomogenousMatrix( float d00, float d01, float d02,
                                float d10, float d11, float d12 )
{
    matrix3_t retval;
    retval.d[0] = d00;
    retval.d[1] = d01;
    retval.d[2] = d02;
    retval.d[3] = d10;
    retval.d[4] = d11;
    retval.d[5] = d12;
    retval.d[6] = 0.0F;
    retval.d[7] = 0.0F;
    retval.d[8] = 1.0F;
    return retval;
}

edge_t MakeEdge( int a, int b )
{
    edge_t retval;
    retval.a = a;
    retval.b = b;
    return retval;
}

vector2_t MakeVector( float x, float y )
{
    vector2_t retval;
    retval.x = x;
    retval.y = y;
    return retval;
}

void Check_vec()
{
    CheckEquiv( normalize3( MakeVector(1,0,0) ), MakeVector(1,0,0) );
    CheckEquiv( normalize3( MakeVector(0,1,0) ), MakeVector(0,1,0) );
    CheckEquiv( normalize3( MakeVector(0,0,1) ), MakeVector(0,0,1) );
    CheckEquiv( normalize3( MakeVector(-1,0,0) ), MakeVector(-1,0,0) );
    CheckEquiv( normalize3( MakeVector(0,-1,0) ), MakeVector(0,-1,0) );
    CheckEquiv( normalize3( MakeVector(0,0,-1) ), MakeVector(0,0,-1) );
    Check( std::abs( norm3( normalize3(MakeVector(1,2,3)) ) - 1.0 ) < 1e-6 );

    CheckEquiv( add3( MakeVector(0,0,0), MakeVector( 4.47F, 3.8F, -6.88F ) ),
                MakeVector( 4.47F, 3.8F, -6.88F ) );
    CheckEquiv( add3( MakeVector( 3, 2, 1 ), MakeVector( -3, -2, -1 ) ),
                MakeVector( 0, 0, 0 ) );
    CheckEquiv( add3( MakeVector( 1, 2, 3 ), MakeVector( -3, -2, -1 ) ),
                MakeVector( -2, 0, 2 ) );

    CheckEquiv( minus3( MakeVector(0,0,0), MakeVector( 4.47F, 3.8F, -6.88F ) ),
                MakeVector( -4.47F, -3.8F, 6.88F ) );
    CheckEquiv( minus3( MakeVector( 3, 2, 1 ), MakeVector( 3, 2, 1 ) ),
                MakeVector( 0, 0, 0 ) );
    CheckEquiv( minus3( MakeVector( 1, 2, 3 ), MakeVector( -3, -2, -1 ) ),
                MakeVector( 4, 4, 4 ) );

    CheckEquiv( scale3( MakeVector(-4.3F,13.9F,4.8e12F ), 0.0F ),
                MakeVector( 0, 0, 0 ) );
    CheckEquiv( scale3( MakeVector(1,2,3), 2 ), MakeVector(2,4,6) );
    CheckEquiv( scale3( MakeVector(1,2,3), -1 ), MakeVector(-1,-2,-3) );

    CheckEquiv( cross3( MakeVector(1,0,0), MakeVector(0,1,0) ),
                MakeVector(0,0,1) );
    CheckEquiv( cross3( MakeVector(0,1,0), MakeVector(1,0,0) ),
                MakeVector(0,0,-1) );
    CheckEquiv( cross3( MakeVector(1,0,0), MakeVector(0,0,1) ),
                MakeVector(0,-1,0) );
    CheckEquiv( cross3( MakeVector(0,0,1), MakeVector(1,0,0) ),
                MakeVector(0,1,0) );
    CheckEquiv( cross3( MakeVector(0,1,0), MakeVector(0,0,1) ),
                MakeVector(1,0,0) );
    CheckEquiv( cross3( MakeVector(0,0,1), MakeVector(0,1,0) ),
                MakeVector(-1,0,0) );
    CheckEquiv( cross3( MakeVector(-1,0,0), MakeVector(0,-1,0) ),
                MakeVector(0,0,1) );
    CheckEquiv( cross3( MakeVector( 3.2F, -4.8F, 0.001F ),
                        MakeVector( 3.2F, -4.8F, 0.001F ) ), MakeVector(0,0,0));

    CheckEquiv( intersect_plane( MakeVector( 1.0F, 0.0F, 0.0F ),
                                 MakeVector( 4.0F, 4.0F, 5.0F ),
                                 MakeVector(-2.0F, 0.0F, 0.0F ) ),
                MakeVector( 0.0F, 4.0F, 5.0F ) );
    CheckEquiv( intersect_plane( MakeVector( 0.0F, 1.0F, 0.0F ),
                                 MakeVector( 4.0F, 4.0F, 5.0F ),
                                 MakeVector( 0.0F, 1.0F, 0.0F ) ),
                MakeVector( 4.0F, 0.0F, 5.0F ) );
    CheckEquiv( intersect_plane( MakeVector( 0.0F, 0.0F, 1.0F ),
                                 MakeVector( 4.0F, 4.0F, 5.0F ),
                                 MakeVector( 0.0F, 0.0F, 4.0F ) ),
                MakeVector( 4.0F, 4.0F, 0.0F ) );

    CheckEquiv( dot3( MakeVector(-3.57F,4.2e-9F,3.2e9F), MakeVector(0,0,0) ),
                0.0F );
    CheckEquiv( dot3( MakeVector(-3.57F,4.2e-9F,0.0), MakeVector(0,0,3.2e9F) ),
                0.0F );
    CheckEquiv( dot3( MakeVector(1,2,3), MakeVector(-3,9,-5) ), 0.0F );
    CheckEquiv( dot3( MakeVector(1,2,3), MakeVector(1,2,3) ), 14.0F );

    CheckEquiv( norm3( MakeVector(0,0,0) ), 0.0F );
    CheckEquiv( norm3( MakeVector(1,0,0) ), 1.0F );
    CheckEquiv( norm3( MakeVector(-1,0,0) ), 1.0F );
    CheckEquiv( norm3( MakeVector(0,1,0) ), 1.0F );
    CheckEquiv( norm3( MakeVector(0,-1,0) ), 1.0F );
    CheckEquiv( norm3( MakeVector(0,0,1) ), 1.0F );
    CheckEquiv( norm3( MakeVector(0,0,-1) ), 1.0F );

    CheckEquiv( multiply3(MakeHomogenousMatrix(1,0,0,0,1,0), 
                          MakeVector(6.4F,-3.9F), 1 ), MakeVector(6.4F,-3.9F));
    CheckEquiv( multiply3(MakeHomogenousMatrix(1,0,3.45F,0,1,-12.62F), 
                          MakeVector(6.4F,-3.9F), 0 ), MakeVector(6.4F,-3.9F));
    CheckEquiv( multiply3(MakeHomogenousMatrix(0.1F,0.2F,0.3F,0.4F,0.5F,0.6F), 
                          MakeVector(1,2), 1 ), MakeVector(0.8F,2.0F), 1E-10F);
    CheckEquiv( multiply3(MakeHomogenousMatrix(0.1F,0.2F,0.3F,0.4F,0.5F,0.6F), 
                          MakeVector(1,2), 0 ), MakeVector(0.5F,1.4F), 1E-10F);

    CheckEquiv( normalize2( MakeVector(1,0) ), MakeVector(1,0) );
    CheckEquiv( normalize2( MakeVector(0,1) ), MakeVector(0,1) );
    CheckEquiv( normalize2( MakeVector(-1,0) ), MakeVector(-1,0) );
    CheckEquiv( normalize2( MakeVector(0,-1) ), MakeVector(0,-1) );
    Check( std::abs( norm2( normalize2(MakeVector(1,2)) ) - 1.0 ) < 1e-6 );

    CheckEquiv( add2( MakeVector(0,0), MakeVector( 4.47F, 3.8F ) ),
                MakeVector( 4.47F, 3.8F ) );
    CheckEquiv( add2( MakeVector( 3, 2 ), MakeVector( -3, -2 ) ),
                MakeVector( 0, 0 ) );
    CheckEquiv( add2( MakeVector( 1, 2 ), MakeVector( -3, -2 ) ),
                MakeVector( -2, 0 ) );

    CheckEquiv( minus2( MakeVector(0,0), MakeVector( 4.47F, 3.8F ) ),
                MakeVector( -4.47F, -3.8F ) );
    CheckEquiv( minus2( MakeVector( 3, 2 ), MakeVector( 3, 2 ) ),
                MakeVector( 0, 0 ) );
    CheckEquiv( minus2( MakeVector( 1, 2 ), MakeVector( -3, -2 ) ),
                MakeVector( 4, 4 ) );

    CheckEquiv( scale2( MakeVector(-4.3F,13.9F ), 0.0F ),
                MakeVector( 0, 0 ) );
    CheckEquiv( scale2( MakeVector(1,2), 2 ), MakeVector(2,4) );
    CheckEquiv( scale2( MakeVector(1,2), -1 ), MakeVector(-1,-2) );

    CheckEquiv( cross2( MakeVector(1,0), MakeVector(0,1) ), 1.0F );
    CheckEquiv( cross2( MakeVector(-1,0), MakeVector(0,1) ), -1.0F );
    CheckEquiv( cross2( MakeVector(1,0), MakeVector(0,-1) ), -1.0F );
    CheckEquiv( cross2( MakeVector(-1,0), MakeVector(0,-1) ), 1.0F );
    CheckEquiv( cross2( MakeVector(0,1), MakeVector(1,0) ), -1.0F );
    CheckEquiv( cross2( MakeVector(0,-1), MakeVector(1,0) ), 1.0F );
    CheckEquiv( cross2( MakeVector(0,1), MakeVector(-1,0) ), 1.0F );
    CheckEquiv( cross2( MakeVector(0,-1), MakeVector(-1,0) ), -1.0F );
    CheckEquiv( cross2( MakeVector( 3.2F, -0.001F ),
                        MakeVector( 3.2F, -0.001F ) ), 0.0F );

    CheckEquiv( dot2( MakeVector(-3.57F,4.2e-9F), MakeVector(0,0) ), 0.0F );
    CheckEquiv( dot2( MakeVector(-3.57F,0.0), MakeVector(0,3.2e9F) ), 0.0F );
    CheckEquiv( dot2( MakeVector(1,2), MakeVector(-4,2) ), 0.0F );
    CheckEquiv( dot2( MakeVector(1,2), MakeVector(1,2) ), 5.0F );

    CheckEquiv( norm2( MakeVector(0,0) ), 0.0F );
    CheckEquiv( norm2( MakeVector(1,0) ), 1.0F );
    CheckEquiv( norm2( MakeVector(-1,0) ), 1.0F );
    CheckEquiv( norm2( MakeVector(0,1) ), 1.0F );
    CheckEquiv( norm2( MakeVector(0,-1) ), 1.0F );

    bool err = 0;
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(-1,0.5F),
                                         MakeVector(1,0), &err), 0.5F, 1E-10 );
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(1,0),
                                         MakeVector(-1,0.5F),&err),0.5F,1E-10);
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(1,0),
                                         MakeVector(-FLT_MIN,1),&err),1.0F,
                1E-10 );
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(-1,0),
                                         MakeVector(-FLT_MIN,1),&err),0.0F,
                1E-10 );
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(1,0),
                                         MakeVector(FLT_MIN,1),&err),0.0F,
                1E-10 );
    CheckEquiv( compute_intersect_scale( MakeVector(0,1), MakeVector(-1,0),
                                         MakeVector(FLT_MIN,1),&err),1.0F,
                1E-10 );
    Check( !err );
    fedisableexcept( FE_UNDERFLOW );
    CheckEquiv( fegetexcept(), FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
    compute_intersect_scale( MakeVector(0,FLT_MIN), MakeVector(-1,0),
                             MakeVector(FLT_MIN,1), &err );
    cout << "DIVBYZERO = " << FE_DIVBYZERO << endl
         << "INVALID = " << FE_INVALID << endl
         << "OVERFLOW = " << FE_OVERFLOW << endl
         << "UNDERFLOW = " << FE_UNDERFLOW << endl
         << "INEXACT = " << FE_INEXACT << endl;
    cout << fetestexcept( FE_ALL_EXCEPT ) << endl;
    feclearexcept( FE_UNDERFLOW );
    feenableexcept( FE_UNDERFLOW );
    cout << fetestexcept( FE_ALL_EXCEPT ) << endl;
    CheckEquiv( fegetexcept(), FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW |
                FE_UNDERFLOW );
    Check( err );
}

void CheckCrossSection
(
 convexpolyhedron_t* poly,
 convexpolygon_t* cs
)
{
    bool* found = new bool[cs->mPointsCnt];
    for( int i = 0; i < cs->mPointsCnt; ++i ) found[i] = false;

    for( int i = 0; i < poly->mPointsCnt; ++i )
    {
        if( std::abs( poly->mPoints[i].z ) < 1e-5F )
        {
            float minDist;
            int minIndex = -1;
            
            for( int j = 0; j < cs->mPointsCnt; ++j )
            {
                float xdist = poly->mPoints[i].x - cs->mPoints[j].x;
                float ydist = poly->mPoints[i].y - cs->mPoints[j].y;
                float dist = sqrt( xdist*xdist + ydist*ydist );
                if( -1 == minIndex || dist < minDist )
                {
                    minIndex = j;
                    minDist = dist;
                }
            }

            Check( -1 != minIndex && std::abs( minDist ) < 1e-5F );
            Check( 0 == found[minIndex] );
            found[minIndex] = 1;
        }
    }

    /* Note: not all cross-section points are found since some are midpoints.*/
    delete[] found;
}

void Check_cut()
{
    dscan_roundcut_measurements_t m;
    m.table_pct = 0.57F;
    m.crownheight_pct = 0.15F;
    m.pavilionheight_pct = 0.445F;
/*    m.girdle_pct = 0.01769F; */
    m.girdle_pct = 0.035F;
    m.star_pct = 0.55F;
    m.pavilion_pct = 0.75F;
    for( int i = 2; i <= 6; ++i )
    { 
        m.girdleedges = i;
        for( int j = 0; j < 2; ++j )
        {
            m.culet_pct = (j == 0 ? 0.0F : 0.015F);
            dscan_roundcut_errors_t errs;
            convexpolyhedron_t* poly = dscan_roundcut_create( m, &errs );
            Check( poly != 0 );
            convexpolyhedron_check( poly, true );
            Check( dscan_roundcut_rotations() > 0 );
            for( int k = 0; k < dscan_roundcut_rotations(); ++k )
            {
                matrix3_t rot = dscan_roundcut_rotation(k);
                convexpolygon_t* cs;
                dscan_convexpolyhedron_transform( poly, rot );
                convexpolyhedron_check( poly, true );
                cs = dscan_roundcut_crosssection(poly,k);
                convexpolygon_check( cs );
                if( 0 == k )
                {
                    CheckCrossSection( poly, cs );
                }
                dscan_convexpolygon_destroy( cs );
            }
            dscan_convexpolyhedron_destroy( poly );
        }
    }

    m.girdle_pct = 0.005F;
    dscan_roundcut_errors_t errs;
    CheckEquiv( dscan_roundcut_create( m, &errs ), (convexpolyhedron_t*)0 );
    Check( errs.min_girdle_pct > m.girdle_pct );
}

void CheckExitCalled( convexpolygon_t poly )
{
    gExitCalled = 0;
    convexpolygon_check( &poly );
    Check( gExitCalled );
    convexpolygon_destroy( poly );
}

void Check_convex_check()
{
    vector2_t* points;

    points = new vector2_t[2];
    points[0].x = points[0].y = 0.0F;
    points[1].x = points[1].y = 1.0F;
    CheckExitCalled( convexpolygon_init( points, 2 ) );

    points = new vector2_t[4];
    points[0] = MakeVector( 0, 0 );
    points[1] = MakeVector( 1, 1 );
    points[2] = MakeVector( 2, 0 );
    points[3] = MakeVector( 2, 1 );
    CheckExitCalled( convexpolygon_init( points, 4 ) );

    points = new vector2_t[6];
    points[0] = MakeVector( 0, 0 );
    points[1] = MakeVector( 2, 0 );
    points[2] = MakeVector( 1, 2 );
    points[3] = MakeVector( 1, 1 );
    points[4] = MakeVector( 2, 2 );
    points[5] = MakeVector( 0, 2 );
    CheckExitCalled( convexpolygon_init( points, 6 ) );
}

convexpolygon_t MakePolygon
(
 float* ps,
 int cnt,
 bool invert
)
{
    if( 0 == cnt ) return convexpolygon_init( 0, 0 );

    vector2_t* vs = (vector2_t*)malloc(sizeof(vector2_t)*cnt);
    if( !invert )
    {
        for( int i = 0; i < cnt; ++i )
        {
            vs[i].x = ps[2*i];
            vs[i].y = ps[2*i+1];
        }
    }
    else
    {
        for( int i = 0; i < cnt; ++i )
        {
            vs[cnt-i-1].x = ps[2*i+1];
            vs[cnt-i-1].y = ps[2*i];
        }
    }
    convexpolygon_t retval = convexpolygon_init( vs, cnt );
    convexpolygon_check( &retval );
    return retval;
}

template<typename T>
convexpolygon_t MakePolygon
(
 const T& first,
 const T& last
)
{
    if( first == last ) return convexpolygon_init( 0, 0 );
    int cnt = std::distance( first, last );
    vector2_t* vs = (vector2_t*)malloc(sizeof(vector2_t)*cnt);

    int i;
    T it;
    for( it = first, i = 0; it != last; ++it, ++i )
    {
        vs[i] = *it;
    }

    convexpolygon_t retval = convexpolygon_init( vs, cnt );
    return retval;
}

float RandomNormal()
{
    float r1 = drand48();
    float r2 = drand48();
    return sqrt( -2 * log(r1) ) * sin( 2 * M_PI * r2 );
}

int RandomPointCnt()
{
    const int MAX = 1000;
    const float MAXLOG = log(MAX);
    const float LOG3 = log(3.0);
    return std::exp( min( MAXLOG, LOG3 + std::abs( RandomNormal() ) ) );
}

vector2_t FromPolar( float r, float theta )
{
    vector2_t retval;
    retval.x = r * cos(theta);
    retval.y = r * sin(theta);
    return retval;
}

struct SwapVec : vector2_t
{
    int index;
    
    SwapVec()
    {
        index = -1;
        x = 0;
        y = 0;
    }

    SwapVec( int i, const vector2_t& impl )
    {
        index = i;
        x = impl.x;
        y = impl.y;
    }
};

struct AngleCompare
{
    vector2_t ref;
    
    bool operator()( const SwapVec& s, const SwapVec& t )
    {
        assert( s.index >= 0 && t.index >= 0 && 
                s.index != t.index );
        if( s.x == ref.x && s.y == ref.y &&
            !(t.x == ref.x && t.y == ref.y ) )
        {
            return 1;
        }
        else if( s.index > t.index )
        {
            return (t.x-ref.x) * (s.y-ref.y) - (t.y-ref.y) * (s.x-ref.x) < 0;
        }
        else
        {
            return (s.x-ref.x) * (t.y-ref.y) - (s.y-ref.y) * (t.x-ref.x) > 0;
        }
    }

    AngleCompare( const vector2_t& s )
    : ref( s )
    {
    }
};

float CounterClockWise
(
 const vector2_t& p1,
 const vector2_t& p2,
 const vector2_t& p3
)
{
    return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

struct VectorCompare
{
    bool operator()( const vector2_t& s, const vector2_t& t )
    {
        if( s.x < t.x ) return 1;
        else if( s.x > t.x ) return 0;
        else return s.y < t.y;
    }
};

vector<vector2_t> RemoveDuplicates( const vector<vector2_t>& in )
{
    vector<vector2_t> retval( in.begin(), in.end() );
    VectorCompare cmp;
    std::sort( retval.begin(), retval.end(), cmp );

    vector<vector2_t>::iterator it = retval.begin();
    while( it + 1 != retval.end() )
    {
        vector<vector2_t>::iterator next = it+1;
        if( it->x == next->x && it->y == next->y )
        {
            it = retval.erase( next ) - 1;
        }
        else
        {
            it = next;
        }
    }

    return retval;
}

convexpolygon_t MakeConvexPolygon
(
 const vector<vector2_t>& dups,
 bool test = 1
)
{
    vector<vector2_t> in = RemoveDuplicates( dups );

    //Graham scan
    vector<SwapVec> s;
    s.reserve( in.size() + 1 );
    vector<vector2_t>::const_iterator it;
    int index;
    for( it = in.begin(), index = 0; it!=in.end(); ++it, ++index )
    {
        s.push_back( SwapVec( index, *it ) );
    }

    int minIndex = 0;
    for( vector<vector2_t>::size_type i = 0; i < s.size(); ++i )
    {
        if( s[i].y < s[minIndex].y ||
            (s[i].y == s[minIndex].y && s[i].x < s[minIndex].x) )
        {
            minIndex = i;
        }
    }

    std::swap( s[minIndex], s[1] );
    s.push_back( s[0] );
    AngleCompare cmp( s[1] );
    std::sort( s.begin()+2, s.end(), cmp );
    s[0] = s[in.size()];

#ifndef NDEBUG
    for( unsigned int i = 1; i < in.size(); ++i )
    {
        assert( !cmp( s[i+1], s[i] ) );
    }
#endif

    int M = 1;
    for( vector<vector2_t>::size_type i = 2; i < s.size(); ++i )
    {
        while( CounterClockWise( s[M-1], s[M], s[i] ) <= 0 )
        {
            if( M > 1 )
            {
                --M;
            }
            else if( i == in.size() )
            {
                break;
            }
            else
            {
                ++i;
            }
        }

        ++M;
        std::swap( s[M], s[i] );
    }

    if( M < 3 ) return convexpolygon_init( 0, 0 );
    if( CounterClockWise( s[M-1], s[0], s[1] ) == 0 ) --M;

    vector2_t* vs = (vector2_t*)malloc(sizeof(vector2_t)*M);
    for( int i = 0; i < M; ++i )
    {
        vs[i] = s[i];
    }
    convexpolygon_t retval = convexpolygon_init( vs, M );
    if( test ) convexpolygon_check( &retval );
    return retval;
}

vector2_t RandomPoint()
{
    vector2_t retval;
    retval.x = drand48()*10 - 5;
    retval.y = drand48()*10 - 5;
    return retval;
}

convexpolygon_t RandomPolygon( bool integerBound )
{

    int cnt = RandomPointCnt();
    vector<vector2_t> points;
    points.reserve(cnt);

    vector2_t offset = FromPolar( drand48(), drand48() * 2.0 * M_PI );
    for( int i = 0; i < cnt; ++i )
    {
        vector2_t point = add2( FromPolar( 1+RandomNormal()*(integerBound?5:1), 
                                           drand48() * 2.0 * M_PI ), offset );
        if( integerBound )
        {
            point.x = floorf( point.x );
            point.y = floorf( point.y );
        }
        points.push_back( point );
    }

    return MakeConvexPolygon( points );
}

matrix3_t GetIdentity()
{
    matrix3_t retval;
    retval.d[0] = 1.0F; retval.d[3] = 0.0F; retval.d[6] = 0.0F;
    retval.d[1] = 0.0F; retval.d[4] = 1.0F; retval.d[7] = 0.0F;
    retval.d[2] = 0.0F; retval.d[5] = 0.0F; retval.d[8] = 1.0F;
    return retval;
}

matrix3_t RandomTransform()
{
    float theta = drand48() * 2 * M_PI;
    float sinTheta = sin(theta);
    float cosTheta = cos(theta);

    matrix3_t rot;
    rot.d[0] = cosTheta; rot.d[1] = -sinTheta; rot.d[2] = drand48() - 0.5F;
    rot.d[3] = sinTheta; rot.d[4] = cosTheta;  rot.d[5] = drand48() - 0.5F;
    rot.d[6] = 0.0F;     rot.d[7] = 0.0F;      rot.d[8] = 1.0F;

    return rot;
}

convexpolygon_t TransformAndCycle
(
 const convexpolygon_t& s, 
 const matrix3_t& rot,
 int* cycle_retval = NULL
)
{
    int cnt = dscan_convexpolygon_points( const_cast<convexpolygon_t*>(&s) );
    if( cnt <= 0 ) return convexpolygon_init( 0, 0 );
    float* rawPoints = new float[cnt*2];
    dscan_convexpolygon_marshall( const_cast<convexpolygon_t*>(&s), rawPoints );
    vector2_t* points = (vector2_t*)malloc(sizeof(vector2_t)*cnt);
    int cycle = (int)(drand48() * cnt);
    assert( cycle < cnt );
    for( int i = 0; i < cnt; ++i )
    {
        float x = rawPoints[2*i];
        float y = rawPoints[2*i + 1];
        points[(i+cycle)%cnt].x = rot.d[0] * x + rot.d[1] * y + rot.d[2];
        points[(i+cycle)%cnt].y = rot.d[3] * x + rot.d[4] * y + rot.d[5];
    }

    convexpolygon_t retval = convexpolygon_init( points, cnt );

    if( cycle_retval != NULL )
    {
        *cycle_retval = cycle;
    }

    delete[] rawPoints;
    return retval;
}

vector<vector2_t> GetPoints( const convexpolygon_t& s )
{
    int cnt = dscan_convexpolygon_points( const_cast<convexpolygon_t*>(&s) );
    if( 0 == cnt ) return vector<vector2_t>();
    float* rawPoints = new float[cnt*2];
    dscan_convexpolygon_marshall( const_cast<convexpolygon_t*>(&s), rawPoints );
    vector<vector2_t> retval;
    retval.reserve( cnt );
    for( int i = 0; i < cnt; ++i )
    {
        vector2_t point;
        point.x = rawPoints[2*i];
        point.y = rawPoints[2*i+1];
        retval.push_back( point );
    }
    delete[] rawPoints;
    return retval;
}

float PointToLine( const vector2_t& point, 
                   const vector2_t& lineStart, const vector2_t& lineEnd )
{
    vector2_t surfaceVec;
    surfaceVec.x = lineEnd.x - lineStart.x;
    surfaceVec.y = lineEnd.y - lineStart.y;

    vector2_t rel;
    rel.x = point.x - lineStart.x;
    rel.y = point.y - lineStart.y;

    float s = norm2( surfaceVec );
    float X = cross2( surfaceVec, rel );
    return X / s;
}

bool PointToLineNear( const vector2_t& point, 
                      const vector2_t& lineStart, const vector2_t& lineEnd,
                      float threshold )
{
    vector2_t surfaceVec;
    surfaceVec.x = lineEnd.x - lineStart.x;
    surfaceVec.y = lineEnd.y - lineStart.y;

    vector2_t rel;
    rel.x = point.x - lineStart.x;
    rel.y = point.y - lineStart.y;

    float s = norm2( surfaceVec );
    float X = cross2( surfaceVec, rel );
    return std::abs(X) < s * threshold;
}

convexpolygon_t Reduce( convexpolygon_t poly, float threshold )
{
    vector<vector2_t> vpoints = GetPoints( poly );
    list<vector2_t> points( vpoints.begin(), vpoints.end() );

    list<vector2_t>::iterator point = points.begin();
    while( point != points.end() && points.size() >= 3 )
    {
        list<vector2_t>::iterator prev;
        if( point != points.begin() ) prev = point; else prev = points.end();
        --prev;
        list<vector2_t>::iterator next = point;
        ++next;
        if( next == points.end() ) next = points.begin();
        float pointDistance = norm2( minus2( *point, *prev ) );
        bool near = PointToLineNear( *point, *prev, *next, threshold );
        if( pointDistance < threshold || near )
        {
            point = points.erase( point );
        }
        else
        {
            ++point;
        }
    }

    convexpolygon_t retval;
    if( points.size() >= 3 )
    {
        
        retval = MakePolygon( points.begin(), points.end() );
    }
    else
    {
        retval = convexpolygon_init( 0, 0 );
    }
    convexpolygon_destroy( poly );

    return retval;
}

void CheckSubset( const convexpolygon_t& s, const convexpolygon_t& t,
                  float threshold )
{
    vector<vector2_t> sp = GetPoints(s);
    vector<vector2_t> tp = GetPoints(t);
    for( vector<vector2_t>::size_type i = 0; i < sp.size(); ++i )
    {
        for( vector<vector2_t>::size_type j = 0; j < tp.size(); ++j )
        {
            Check( PointToLine( sp[i], tp[j], tp[(j+1)%tp.size()] ) >
                   -threshold );
        }
    }
}

void CheckIntersect( const convexpolygon_t& s, const convexpolygon_t& t,
                     const convexpolygon_t& result )
{
    for( int i = 0; i < XITERS; ++i )
    {
        matrix3_t mat = (0 == i ? GetIdentity() : RandomTransform() );
        convexpolygon_t sp = TransformAndCycle( s, mat );
        convexpolygon_t tp = TransformAndCycle( t, mat );
        convexpolygon_t x = Reduce( convexpolygon_intersect(
                                     &const_cast<convexpolygon_t&>(sp),
                                     &const_cast<convexpolygon_t&>(tp) ),
                                    1e-4 );
        convexpolygon_t test = TransformAndCycle( result, mat );
        convexpolygon_check( &x );
        convexpolygon_check( &test );
        CheckEquiv( x, test, 1e-3F );
        convexpolygon_destroy( sp );
        convexpolygon_destroy( tp );
        convexpolygon_destroy( x );
        convexpolygon_destroy( test );
    }
}

void CheckIntersect( const convexpolygon_t& s, const convexpolygon_t& t )
{
    const float THRESHOLD = 1e-4;
    convexpolygon_t raw = convexpolygon_intersect(
     &const_cast<convexpolygon_t&>(s), &const_cast<convexpolygon_t&>(t) );
    CheckSubset( raw, s, THRESHOLD );
    CheckSubset( raw, t, THRESHOLD );
    convexpolygon_t x = Reduce( raw, THRESHOLD );
    CheckIntersect( s, t, x );
    convexpolygon_destroy( x );
}

void CheckIntersect
(
 float ps[],
 int ls,
 float pt[],
 int lt,
 float pr[],
 int lr
)
{
    for( int i = 0; i < 2; ++i )
    {
        convexpolygon_t s = MakePolygon( ps, ls, (i==1) );
        convexpolygon_t t = MakePolygon( pt, lt, (i==1) );
        convexpolygon_t r = MakePolygon( pr, lr, (i==1) );
        CheckIntersect( s, t, r );
        CheckIntersect( t, s, r );
        convexpolygon_destroy( s );
        convexpolygon_destroy( t );
        convexpolygon_destroy( r );
    }
}

void Check_convex_intersect()
{
    /* No intersection cases */
    // No comparison
    float ps0[] = {0,0, 1,0, 1,1, 0,1};
    float* pt0 = NULL;
    float* pr0 = NULL;
    CheckIntersect( ps0, 4, pt0, 0, pr0, 0 );
    CheckIntersect( pt0, 0, ps0, 4, pr0, 0 );
    CheckIntersect( pt0, 0, pt0, 0, pr0, 0 );

    // No parallels
    float ps1[] = {-1,2, -1,1, 0,1, 0,2};
    float pt1[] = {1,-1, 1,-2, 2,-2, 2,-1};
    float* pr1 = NULL;
    CheckIntersect( ps1, 4, pt1, 4, pr1, 0 );

    // Directed parallels
    float ps2[] = {-2,-1, -1,-1, -1, 1, -2, 1};
    float pt2[] = {2,1, 1,1, 1,-1, 2, -1};
    float* pr2 = NULL;
    CheckIntersect( ps2, 4, pt2, 4, pr2, 0 );

    // Opposing parallels
    float ps3[] = {-1,2, -1,1, 0,1, 0,2};
    float pt3[] = {0,-1, 0,-2, 1,-2, 1,-1};
    float* pr3 = NULL;
    CheckIntersect( ps3, 4, pt3, 4, pr3, 0 );

    /* Degenerate intersection cases */
    // Point to point
    float ps4[] = {0,1, 1,1, 1,2, 0,2};
    float pt4[] = {1,0, 2,0, 2,1, 1,1};
    float* pr4 = NULL;
    CheckIntersect( ps4, 4, pt4, 4, pr4, 0 );

    // Point to edge
    float ps5[] = {0,2, 2,2, 2,4, 0,4};
    float pt5[] = {1,0, 2,1, 1,2, 0,1};
    float* pr5 = NULL;
    CheckIntersect( ps5, 4, pt5, 4, pr5, 0 );

    // Edge to edge and one point to edge
    float ps6[] = {1,0, 3,0, 3,1, 1,1};
    float pt6[] = {0,1, 2,1, 2,2, 0,2};
    float* pr6 = NULL;
    CheckIntersect( ps6, 4, pt6, 4, pr6, 0 );

    // Edge to edge and two point to edge
    float ps7[] = {1,0, 3,0, 3,1, 1,1};
    float pt7[] = {0,1, 4,1, 4,2, 0,2};
    float* pr7 = NULL;
    CheckIntersect( ps7, 4, pt7, 4, pr7, 0 );

    // Edge to edge and point to point
    float ps8[] = {1,0, 3,0, 3,1, 1,1};
    float pt8[] = {1,1, 3,1, 3,2, 1,2};
    float* pr8 = NULL;
    CheckIntersect( ps8, 4, pt8, 4, pr8, 0 );

    /* Intersection cases */
    // Contained
    float ps9[] = {2,2, 3,2, 3,3, 2,3};
    float pt9[] = {1,1, 4,1, 4,4, 1,4};
    CheckIntersect( ps9, 4, pt9, 4, ps9, 4 );

    // Edge to edge
    float ps10[] = {25,25, 50,50, 25,75, 0,50};
    float pt10[] = {49,25, 74,50, 49,75, 24,50};
    float pr10[] = {37,63, 24,50, 37,37, 50,50};
    CheckIntersect( ps10, 4, pt10, 4, pr10, 4 );

    // Dual edge to edge
    float ps11[] = {0,2, 2,0, 4,2, 2,4};
    float pt11[] = {3,0, 7,0, 7,4, 3,4};
    float pr11[] = {3,1, 4,2, 3,3};
    CheckIntersect( ps11, 4, pt11, 4, pr11, 3 );

    // Point to edge
    float ps12[] = {0,1, -1,0, 0,-1, 1,0};
    float pt12[] = {2,2, -2,2, -2,0, 2,0};
    float pr12[] = {0,1, -1,0, 1,0};
    CheckIntersect( ps12, 4, pt12, 4, pr12, 3 );

    // Point to point
    float ps13[] = {-1,1, -2,0, -2,-2, -1,-3, 1,-3, 2,-2, 2,0, 1,1};
    float pt13[] = {1,-1, 2,0, 2,2, 1,3, -1,3, -2,2, -2,0, -1,-1};
    float pr13[] = {-1,1, -2,0, -1,-1, 1,-1, 2,0, 1,1};
    CheckIntersect( ps13, 8, pt13, 8, pr13, 6 );

    // Point to point and point to edge
    float ps14[] = {0,5, -5,0, 0,-5, 5,0};
    float pt14[] = {0,5, 0,-5, 5,-10, 5,10};
    float pr14[] = {0,5, 0,-5, 5,0};
    CheckIntersect( ps14, 4, pt14, 4, pr14, 3 );

    // Parallel edges
    float ps15[] = {-2,-2, 2,-2, 2,2, -2,2};
    float pt15[] = {-1,-2, 1,-2, 1,2, -1,2};
    CheckIntersect( ps15, 4, pt15, 4, pt15, 4 );

    for( int j = 0; j < 2; ++j )
    {
        for( int i = 0; i < CVXITERS; ++i )
        {
            convexpolygon_t s = RandomPolygon( (j==1) );
            convexpolygon_t t = RandomPolygon( (j==1) );
            CheckIntersect( s, t );
            convexpolygon_destroy( s );
            convexpolygon_destroy( t );
        }
        cout << "OK\n";
    }
}

convexpolygon_t GetProjectIntersect
(
 convexpolygon_t& poly,
 vector2_t surface_start,
 vector2_t surface_end,
 vector2_t ray,
 float tolerance
)
{
    projection_t r0 = convexpolygon_project( &poly, surface_start, surface_end,
                                             ray );
    projection_t r1 = convexpolygon_project( &poly, surface_end, surface_start,
                                             scale2( ray, -1.0F ) );

    if( r0.edge_start >= 0 )
    {
        Check( r1.edge_start >= 0 );
        vector<vector2_t> points;

        points.push_back( add2( poly.mPoints[r0.edge_start],
                                scale2( poly.mVecs[r0.edge_start],
                                        r0.scale_start ) ) );
        for( int i = r0.edge_start; i != r0.edge_end; i = (i+1)%poly.mPointsCnt )
        {
            points.push_back( poly.mPoints[(i+1)%poly.mPointsCnt] );
        }
        points.push_back( add2( poly.mPoints[r0.edge_end],
                                scale2( poly.mVecs[r0.edge_end],
                                        r0.scale_end ) ) );

        points.push_back( add2( poly.mPoints[r1.edge_start],
                                scale2( poly.mVecs[r1.edge_start],
                                        r1.scale_start ) ) );
        for( int i = r1.edge_start; i != r1.edge_end; i = (i+1)%poly.mPointsCnt )
        {
            points.push_back( poly.mPoints[(i+1)%poly.mPointsCnt] );
        }
        points.push_back( add2( poly.mPoints[r1.edge_end],
                                scale2( poly.mVecs[r1.edge_end],
                                        r1.scale_end ) ) );
        return Reduce( MakeConvexPolygon( points, false ), tolerance );
    }
    else
    {
        Check( r1.edge_start < 0 );
        return convexpolygon_init( NULL, 0 );
    }
}

convexpolygon_t GetIntersect
(
 convexpolygon_t& poly,
 vector2_t surface_start,
 vector2_t surface_end,
 vector2_t ray
)
{
    float maxDist = max( norm2( minus2( poly.center, surface_start ) ),
                         norm2( minus2( poly.center, surface_end ) ) );
    float radius = 2 * (maxDist + sqrt(poly.radius2));
    vector<vector2_t> points;
    vector2_t nray = normalize2( ray );
    points.push_back( MakeVector(
                       surface_start.x + scale2( nray, radius ).x,
                       surface_start.y + scale2( nray, radius ).y ) );
    points.push_back( MakeVector( 
                       surface_start.x + scale2( nray, -radius ).x,
                       surface_start.y + scale2( nray, -radius ).y ) );
    points.push_back( MakeVector(
                       surface_end.x + scale2( nray, -radius ).x,
                       surface_end.y + scale2( nray, -radius ).y ) );
    points.push_back( MakeVector( 
                       surface_end.x + scale2( nray, radius ).x,
                       surface_end.y + scale2( nray, radius ).y ) );
    convexpolygon_t area = MakeConvexPolygon( points );
    return convexpolygon_intersect( &poly, &area );
}

static vector2_t HMultiply(matrix3_t m, vector2_t v)
{
    vector2_t retval;
    retval.x = m.d[0] * v.x + m.d[1] * v.y + m.d[2];
    retval.y = m.d[3] * v.x + m.d[4] * v.y + m.d[5];
    return retval;
}

static vector2_t NonHMultiply(matrix3_t m, vector2_t v)
{
    vector2_t retval;
    retval.x = m.d[0] * v.x + m.d[1] * v.y;
    retval.y = m.d[3] * v.x + m.d[4] * v.y;
    return retval;
}

void CheckProject
(
 convexpolygon_t& poly,
 vector2_t surface_start,
 vector2_t surface_end,
 vector2_t ray,
 convexpolygon_t& test
)
{
    for( int i = 0; i < XITERS; ++i )
    {
        matrix3_t mat = (0 == i ? GetIdentity() : RandomTransform() );
        convexpolygon_t polyPrime = TransformAndCycle( poly, mat );
        vector2_t startPrime = HMultiply( mat, surface_start );
        vector2_t endPrime = HMultiply( mat, surface_end );
        vector2_t rayPrime = NonHMultiply( mat, ray );
        convexpolygon_t testPrime = TransformAndCycle( test, mat );
        convexpolygon_t x = GetProjectIntersect( polyPrime, startPrime,
                                                 endPrime, rayPrime,
                                                 1e-4F );
        CheckEquiv( x, testPrime, 1e-4F );
        convexpolygon_destroy( polyPrime );
        convexpolygon_destroy( testPrime );
    }
}

projection_t Cycle( projection_t p, convexpolygon_t& poly, int cycle )
{
    projection_t retval = p;
    if( p.edge_start >= 0 )
    {
        retval.edge_start = (p.edge_start + poly.mPointsCnt + cycle) %
            poly.mPointsCnt;
    }
    if( p.edge_end >= 0 )
    {
        retval.edge_end = (p.edge_end + poly.mPointsCnt + cycle) %
            poly.mPointsCnt;
    }
    return retval;
}

void CheckProject
(
 float data[], 
 int sides,
 vector2_t start, 
 vector2_t end,
 vector2_t ray,
 int edge0,
 float scale0,
 int edge1,
 float scale1,
 float tolerance
)
{
    vector2_t* points = (vector2_t*) malloc( sizeof(vector2_t) * sides );
    for( int i = 0; i < sides; ++i )
    {
        points[i] = MakeVector( data[2*i], data[2*i + 1] );
    }
    convexpolygon_t poly = convexpolygon_init( sides ? points : NULL, sides );
    if( 0 == sides ) free( points );

    for( int i = 0; i < XITERS; ++i )
    {
        matrix3_t mat = (0 == i ? GetIdentity() : RandomTransform() );
        int cycle;
        convexpolygon_t polyPrime = TransformAndCycle( poly, mat, &cycle );
        vector2_t startPrime = HMultiply( mat, start );
        vector2_t endPrime = HMultiply( mat, end );
        vector2_t rayPrime = NonHMultiply( mat, ray );
        projection_t r = Cycle( convexpolygon_project( &polyPrime, startPrime,
                                                       endPrime, rayPrime ),
                                polyPrime, -cycle );
        CheckEquiv( r.edge_start, edge0 );
        CheckEquiv( r.scale_start, scale0, tolerance );
        CheckEquiv( r.edge_end, edge1 );
        CheckEquiv( r.scale_end, scale1, tolerance );
        convexpolygon_destroy( polyPrime );
    }

    convexpolygon_destroy( poly );
}

void CheckProject
(
 convexpolygon_t& poly,
 vector2_t surface_start,
 vector2_t surface_end
)
{
    vector2_t vec = minus2( surface_end, surface_start );
    vector2_t ray;
    ray.x = -vec.y;
    ray.y = vec.x;
    convexpolygon_t test = GetIntersect( poly, surface_start, surface_end, 
                                         ray );
    CheckProject( poly, surface_start, surface_end, ray, test );
    convexpolygon_destroy( test );
}

void CheckProject
(
 float data[],
 int sides,
 vector2_t surface_start,
 vector2_t surface_end
)
{
    vector2_t* points = (vector2_t*) malloc( sizeof(vector2_t) * sides );
    for( int i = 0; i < sides; ++i )
    {
        points[i] = MakeVector( data[2*i], data[2*i + 1] );
    }
    convexpolygon_t poly = convexpolygon_init( points, sides );

    CheckProject( poly, surface_start, surface_end );
    convexpolygon_destroy( poly );
}

void CheckProject
(
 float data[],
 int sides,
 vector2_t surface_start,
 vector2_t surface_end,
 vector2_t ray
)
{
    vector<vector2_t> points;
    for( int i = 0; i < sides; ++i )
    {
        points.push_back( MakeVector( data[2*i], data[2*i + 1] ) );
    }
    convexpolygon_t poly = MakeConvexPolygon( points );

    convexpolygon_t test = GetIntersect( poly, surface_start, surface_end, 
                                         ray );
    CheckProject( poly, surface_start, surface_end, ray, test );
    convexpolygon_destroy( test );
    convexpolygon_destroy( poly );
}

void Check_convex_project()
{
    float* p0 = NULL;
    vector2_t s0 = MakeVector( -2, -1 );
    vector2_t e0 = MakeVector( -1, -1 );
    vector2_t r0 = MakeVector(  0,  1 );
    vector2_t n0 = MakeVector( 0, -1 );
    CheckProject( p0, 0, s0, e0, r0, -1, 0, -1, 0, 0 );
    CheckProject( p0, 0, e0, s0, n0, -1, 0, -1, 0, 0 );

    float p1[] = {0, 0, 2, 0, 1, 1 };
    vector2_t s1 = MakeVector( -2, -1 );
    vector2_t e1 = MakeVector( -1, -1 );
    vector2_t r1 = MakeVector(  0,  1 );
    vector2_t n1 = MakeVector( 0, -1 );
    CheckProject( p1, 3, s1, e1, r1, -1, 0, -1, 0, 0 );
    CheckProject( p1, 3, e1, s1, n1, -1, 0, -1, 0, 0 );

    float p2[] = {0, 0, 2, 0, 1, 1 };
    vector2_t s2 = MakeVector( 3, 3 );
    vector2_t e2 = MakeVector( 4, 3 );
    vector2_t r2 = MakeVector(  0,  1 );
    vector2_t n2 = MakeVector( 0, -1 );
    CheckProject( p2, 3, s2, e2, r2, -1, 0, -1, 0, 0 );
    CheckProject( p2, 3, e2, s2, n2, -1, 0, -1, 0, 0 );

    float p3[] = {0, 0, 2, 0, 1, 1 };
    vector2_t s3 = MakeVector( -2, -1 );
    vector2_t e3 = MakeVector(  0.5F, -1 );
    vector2_t r3 = MakeVector(  0,  1 );
    vector2_t n3 = MakeVector( 0, -1 );
    CheckProject( p3, 3, s3, e3, r3, 0, 0, 0, 0.25F, 1E-05F );
    CheckProject( p3, 3, e3, s3, n3, 2, 0.5F, 2, 1, 1E-05F );

    float p4[] = {0, 0, 2, 0, 1, 1 };
    vector2_t s4 = MakeVector(  1.5F, -1 );
    vector2_t e4 = MakeVector(  3, -1 );
    vector2_t r4 = MakeVector(  0,  1 );
    vector2_t n4 = MakeVector( 0, -1 );
    CheckProject( p4, 3, s4, e4, r4, 0, 0.75F, 0, 1, 1E-05F );
    CheckProject( p4, 3, e4, s4, n4, 1, 0, 1, 0.5F, 1E-05F );

    float p5[] = {0, 0, 2, 0, 1, 1 };
    vector2_t s5 = MakeVector(  0.5F, -1 );
    vector2_t e5 = MakeVector(  1.5F, -1 );
    vector2_t r5 = MakeVector(  0,  1 );
    vector2_t n5 = MakeVector( 0, -1 );
    CheckProject( p5, 3, s5, e5, r5, 0, 0.25F, 0, 0.75F, 1E-05F );
    CheckProject( p5, 3, e5, s5, n5, 1, 0.5F, 2, 0.5F, 1E-05F );

    float p6[] = {0, 0, 1, 0, 1, 1, 0, 1 };
    vector2_t s6 = MakeVector(  0, -1 );
    vector2_t e6 = MakeVector(  1, -1 );
    CheckProject( p6, 4, s6, e6 );

    float p7[] = {0, 0, 1, 0, 1, 1, 0, 1 };
    vector2_t s7 = MakeVector(  -1, -1 );
    vector2_t e7 = MakeVector(  2, -1 );
    CheckProject( p7, 4, s7, e7 );

    float p8[] = {0, 0, 1, 0, 1, 1, 0, 1 };
    vector2_t s8 = MakeVector(  0.25F, -1 );
    vector2_t e8 = MakeVector(  0.75F, -1 );
    CheckProject( p8, 4, s8, e8 );

    float p9[] = {0, 0, 1, 0, 1, 1, 0, 1 };
    vector2_t s9 = MakeVector(  0.5F, -1 );
    vector2_t e9 = MakeVector(  0.5F, -1 );
    vector2_t r9 = MakeVector( 0, 1 );
    CheckProject( p9, 4, s9, e9, r9 );

    float p10[] = {0, 0, 1, 0, 1, 1, 0, 1 };
    vector2_t s10 = MakeVector(  0.0F, -1 );
    vector2_t e10 = MakeVector(  0.0F, -1 );
    vector2_t r10 = MakeVector( 0, 1 );
    CheckProject( p10, 4, s10, e10, r10 );

    for( int j = 0; j < 2; ++j )
    {
        for( int i = 0; i < CVXITERS; ++i )
        {
            convexpolygon_t s = Reduce( RandomPolygon( (j==1) ), 1E-04F );
            vector2_t sg0 = RandomPoint();
            vector2_t sg1 = RandomPoint();
            if( norm2( minus2( sg0, sg1 ) ) > 1E-4F )
            {
                CheckProject( s, sg0, sg1 );
            }
            convexpolygon_destroy( s );
        }
    }
}

static vector3_t MakeVector3( float x, float y, float z )
{
    vector3_t retval;
    retval.x = x;
    retval.y = y;
    retval.z = z;
    return retval;
}

static face_t MakeFace3( int e0, int e1, int e2, bool r0, bool r1 )
{
    face_t retval;
    retval.edges = (int*)malloc(sizeof(int)*3);
    retval.edges[0] = e0;
    retval.edges[1] = e1;
    retval.edges[2] = e2;
    retval.edgesCnt = 3;
    retval.reverse0 = r0;
    retval.reverse1 = r1;
    return retval;
}

static face_t MakeFace4( int e0, int e1, int e2, int e3, bool r0, bool r1 )
{
    face_t retval;
    retval.edges = (int*)malloc(sizeof(int)*4);
    retval.edges[0] = e0;
    retval.edges[1] = e1;
    retval.edges[2] = e2;
    retval.edges[3] = e3;
    retval.edgesCnt = 4;
    retval.reverse0 = r0;
    retval.reverse1 = r1;
    return retval;
}

dscan_convexpolyhedron_t MakeRegularTetrahedron()
{
    convexpolyhedron_t* retval = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));

    vector3_t* points = (vector3_t*)malloc(sizeof(vector3_t)*4);
    edge_t* edges = (edge_t*)malloc(sizeof(edge_t)*6);
    face_t* faces = (face_t*)malloc(sizeof(face_t)*4);

    points[0] = MakeVector3( -0.5F, -sqrt(3)/6.0F,         0.0F );
    points[1] = MakeVector3(  0.5F, -sqrt(3)/6.0F,         0.0F );
    points[2] = MakeVector3(  0.0F,  sqrt(3)/3.0F,         0.0F );
    points[3] = MakeVector3(  0.0F,          0.0F, sqrt(6)/3.0F );

    edges[0].a = 0;
    edges[0].b = 1;
    edges[1].a = 1;
    edges[1].b = 2;
    edges[2].a = 2;
    edges[2].b = 0;
    edges[3].a = 0;
    edges[3].b = 3;
    edges[4].a = 1;
    edges[4].b = 3;
    edges[5].a = 2;
    edges[5].b = 3;

    faces[0] = MakeFace3( 2, 1, 0,  1, 1 );
    faces[1] = MakeFace3( 0, 4, 3,  0, 0 );
    faces[2] = MakeFace3( 1, 5, 4,  0, 0 );
    faces[3] = MakeFace3( 2, 3, 5,  0, 0 );

    *retval = convexpolyhedron_init( points, 4, edges, 6, faces, 4 );
    convexpolyhedron_check( retval, 1 );
    return retval;
}

struct SortableEdge
{
    float x0;
    float y0;
    float x1;
    float y1;

    SortableEdge( float x0i, float y0i, float x1i, float y1i )
    {
        if( x0i < x1i || (x0i == x1i && y0i < y1i) )
        {
            x0 = x0i;
            y0 = y0i;
            x1 = x1i;
            y1 = y1i;
        }
        else
        {
            x1 = x0i;
            y1 = y0i;
            x0 = x1i;
            y0 = y1i;
        }
    }

    bool operator<( const SortableEdge& that ) const
    {
        if( x0 < that.x0 ) return 1;
        else if( x0 > that.x0 ) return 0;
        else if( y0 < that.y0 ) return 1;
        else if( y0 > that.y0 ) return 0;
        else if( x1 < that.x1 ) return 1;
        else if( x1 > that.x1 ) return 0;
        else return y1 < that.y1;
    }
};

SortableEdge MakeSortableEdge( dscan_convexpolyhedron_t poly, int edge )
{
    return SortableEdge( poly->mPoints[poly->mEdges[edge].a].x,
                         poly->mPoints[poly->mEdges[edge].a].y,
                         poly->mPoints[poly->mEdges[edge].b].x,
                         poly->mPoints[poly->mEdges[edge].b].y );                         
}

void CheckTetraVisible( dscan_convexpolyhedron_t tetra, bool v0, bool v1,
                        bool v2, bool v3, bool v4, bool v5 )
{
    float* edgePoints = new float[dscan_convexpolyhedron_edges(tetra)*4];
    int cnt = dscan_convexpolyhedron_visible(tetra, edgePoints);
    set<SortableEdge> expected;
    if( v0 ) expected.insert( MakeSortableEdge( tetra, 0 ) );
    if( v1 ) expected.insert( MakeSortableEdge( tetra, 1 ) );
    if( v2 ) expected.insert( MakeSortableEdge( tetra, 2 ) );
    if( v3 ) expected.insert( MakeSortableEdge( tetra, 3 ) );
    if( v4 ) expected.insert( MakeSortableEdge( tetra, 4 ) );
    if( v5 ) expected.insert( MakeSortableEdge( tetra, 5 ) );
    set<SortableEdge> found;
    for( int i = 0; i < cnt; i += 4 )
    {
        found.insert( SortableEdge( edgePoints[i+0], edgePoints[i+1],
                                    edgePoints[i+2], edgePoints[i+3] ) );
    }
    CheckEquiv( expected, found );
    delete[] edgePoints;
}

matrix3_t RandomTransform3()
{
    float theta = drand48() * 2 * M_PI;
    float phi = drand48() * 2 * M_PI;
    float psi = drand48() * 2 * M_PI;

    matrix3_t rot;
    rot.d[0] = cos(theta)*cos(psi); 
    rot.d[1] = cos(theta)*sin(psi);
    rot.d[2] = -sin(theta);

    rot.d[3] = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
    rot.d[4] = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
    rot.d[5] = sin(phi)*cos(theta);

    rot.d[6] = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
    rot.d[7] = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
    rot.d[8] = cos(phi)*cos(theta);

    return rot;
}

matrix3_t MakeXRotation3( float degrees )
{
    float theta = degrees * M_PI / 180;

    matrix3_t rot;
    rot.d[0] = 1.0; 
    rot.d[1] = 0.0;
    rot.d[2] = 0.0;

    rot.d[3] = 0.0;
    rot.d[4] = cos(theta);
    rot.d[5] = -sin(theta);

    rot.d[6] = 0.0;
    rot.d[7] = sin(theta);
    rot.d[8] = cos(theta);

    return rot;
}

void CheckExitCalled( dscan_convexpolyhedron_t poly )
{
    gExitCalled = 0;
    convexpolyhedron_check( poly, 0 );
    Check( gExitCalled );
    dscan_convexpolyhedron_destroy( poly );
}

void Check_convex_3D_check()
{
    convexpolyhedron_t* poly;
    vector3_t* points;
    edge_t* edges;
    face_t* faces;

    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( NULL, 0, NULL, 0, NULL, 0 );
    CheckExitCalled( poly );

    points = (vector3_t*) malloc( sizeof(vector3_t) * 1 );
    points[0] = MakeVector( 0, 0, 0 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 1, NULL, 0, NULL, 0 );
    CheckExitCalled( poly );

    points = (vector3_t*) malloc( sizeof(vector3_t) * 2 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 1, 1, 1 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 1 );
    edges[0] = MakeEdge( 0, 1 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 2, edges, 1, NULL, 0 );
    CheckExitCalled( poly );

    //Degenerate Face
    points = (vector3_t*) malloc( sizeof(vector3_t) * 2 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 1 );
    faces = (face_t*) malloc( sizeof(face_t) * 1 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 1, 1, 1 );
    edges[0] = MakeEdge( 0, 1 );
    faces[0].edges = (int*)malloc(sizeof(int)*1);
    faces[0].edges[0] = 0;
    faces[0].edgesCnt = 1;
    faces[0].reverse0 = 0;
    faces[0].reverse1 = 0;
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 2, edges, 1, faces, 1 );
    CheckExitCalled( poly );

    //Discontinuous Face
    points = (vector3_t*) malloc( sizeof(vector3_t) * 4 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 3 );
    faces = (face_t*) malloc( sizeof(face_t) * 2 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 0, 1, 0 );
    points[2] = MakeVector( 1, 0, 0 );
    points[3] = MakeVector( 0, 0, 1 );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 2, 1 );
    edges[2] = MakeEdge( 0, 3 );
    faces[0] = MakeFace3( 0, 1, 2, 0, 1 );
    faces[1] = MakeFace3( 2, 1, 0, 0, 0 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 4, edges, 3, faces, 2 );
    CheckExitCalled( poly );

    //Non-colinear face
    points = (vector3_t*) malloc( sizeof(vector3_t) * 5 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 8 );
    faces = (face_t*) malloc( sizeof(face_t) * 5 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 1, 0, 0 );
    points[2] = MakeVector( 1, 1, 0.25F );
    points[3] = MakeVector( 0, 1, 0 );
    points[4] = MakeVector( 0.5F, 0.5F, 1 );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 1, 2 );
    edges[2] = MakeEdge( 2, 3 );
    edges[3] = MakeEdge( 3, 0 );
    edges[4] = MakeEdge( 0, 4 );
    edges[5] = MakeEdge( 1, 4 );
    edges[6] = MakeEdge( 2, 4 );
    edges[7] = MakeEdge( 3, 4 );
    faces[0] = MakeFace4( 0, 1, 2, 3, 0, 0 );
    faces[1] = MakeFace3( 0, 5, 4, 0, 0 );
    faces[2] = MakeFace3( 1, 6, 5, 0, 0 );
    faces[3] = MakeFace3( 2, 7, 6, 0, 0 );
    faces[4] = MakeFace3( 3, 4, 7, 0, 0 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 5, edges, 8, faces, 5 );
    CheckExitCalled( poly );

    //Non-convex
    points = (vector3_t*) malloc( sizeof(vector3_t) * 9 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 16 );
    faces = (face_t*) malloc( sizeof(face_t) * 9 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 1, 0, 0 );
    points[2] = MakeVector( 1, 1, 0 );
    points[3] = MakeVector( 0, 1, 0 );
    points[4] = MakeVector( 0, 0, 1 );
    points[5] = MakeVector( 1, 0, 1 );
    points[6] = MakeVector( 1, 1, 1 );
    points[7] = MakeVector( 0, 1, 1 );
    points[8] = MakeVector( 0.5F, 0.5F, 0.5F );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 1, 2 );
    edges[2] = MakeEdge( 2, 3 );
    edges[3] = MakeEdge( 3, 0 );
    edges[4] = MakeEdge( 0, 4 );
    edges[5] = MakeEdge( 1, 5 );
    edges[6] = MakeEdge( 2, 6 );
    edges[7] = MakeEdge( 3, 7 );
    edges[8] = MakeEdge( 4, 5 );
    edges[9] = MakeEdge( 5, 6 );
    edges[10] = MakeEdge( 6, 7 );
    edges[11] = MakeEdge( 7, 4 );
    edges[12] = MakeEdge( 4, 8 );
    edges[13] = MakeEdge( 5, 8 );
    edges[14] = MakeEdge( 6, 8 );
    edges[15] = MakeEdge( 7, 8 );
    faces[0] = MakeFace4( 3, 2, 1, 0, 1, 1 );
    faces[1] = MakeFace4( 0, 5, 8, 4, 0, 0 );
    faces[2] = MakeFace4( 1, 6, 9, 5, 0, 0 );
    faces[3] = MakeFace4( 2, 7, 10, 6, 0, 0 );
    faces[4] = MakeFace4( 3, 4, 11, 7, 0, 0 );
    faces[5] = MakeFace3( 8, 13, 12, 0, 0 );
    faces[6] = MakeFace3( 9, 14, 13, 0, 0 );
    faces[7] = MakeFace3( 10, 15, 14, 0, 0 );
    faces[8] = MakeFace3( 11, 12, 15, 0, 0 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 9, edges, 16, faces, 9 );
    CheckExitCalled( poly );

    //Extra edge
    points = (vector3_t*) malloc( sizeof(vector3_t) * 3 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 4 );
    faces = (face_t*) malloc( sizeof(face_t) * 2 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 0, 1, 0 );
    points[2] = MakeVector( 1, 0, 0 );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 1, 2 );
    edges[2] = MakeEdge( 2, 0 );
    edges[3] = MakeEdge( 2, 0 );
    faces[0] = MakeFace3( 0, 1, 2, 0, 0 );
    faces[1] = MakeFace3( 2, 1, 0, 1, 1 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 3, edges, 4, faces, 2 );
    CheckExitCalled( poly );

    points = (vector3_t*) malloc( sizeof(vector3_t) * 3 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 4 );
    faces = (face_t*) malloc( sizeof(face_t) * 2 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 0, 1, 0 );
    points[2] = MakeVector( 1, 0, 0 );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 1, 2 );
    edges[2] = MakeEdge( 2, 0 );
    edges[3] = MakeEdge( 2, 0 );
    faces[0] = MakeFace3( 0, 1, 2, 0, 0 );
    faces[1] = MakeFace3( 3, 1, 0, 1, 1 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 3, edges, 4, faces, 2 );
    CheckExitCalled( poly );

    //Extra point
    points = (vector3_t*) malloc( sizeof(vector3_t) * 4 );
    edges = (edge_t*) malloc( sizeof(edge_t) * 3 );
    faces = (face_t*) malloc( sizeof(face_t) * 2 );
    points[0] = MakeVector( 0, 0, 0 );
    points[1] = MakeVector( 0, 1, 0 );
    points[2] = MakeVector( 1, 0, 0 );
    points[3] = MakeVector( 0, 0, 1 );
    edges[0] = MakeEdge( 0, 1 );
    edges[1] = MakeEdge( 1, 2 );
    edges[2] = MakeEdge( 2, 0 );
    faces[0] = MakeFace3( 0, 1, 2, 0, 0 );
    faces[1] = MakeFace3( 2, 1, 0, 1, 1 );
    poly = (convexpolyhedron_t*)malloc(sizeof(convexpolyhedron_t));
    *poly = convexpolyhedron_init( points, 4, edges, 3, faces, 2 );
    CheckExitCalled( poly );
}

void Check_convex_3D()
{
    dscan_convexpolyhedron_t tetra = MakeRegularTetrahedron();
    CheckEquiv( dscan_convexpolyhedron_volume( tetra ), 
                (float)sqrt(2)/12.0F, 1E-4F );
    CheckEquiv( dscan_convexpolyhedron_edges( tetra ), 6 );
    CheckTetraVisible( tetra, 1, 1, 1, 1, 1, 1 );

    dscan_convexpolyhedron_transform( tetra, MakeXRotation3( 54.7356103F ) );
    CheckTetraVisible( tetra, 0, 1, 1, 1, 1, 1 );

    dscan_convexpolyhedron_transform( tetra, RandomTransform3() );
    convexpolyhedron_check( tetra, 1 );
    CheckEquiv( dscan_convexpolyhedron_volume( tetra ), 
                (float)sqrt(2)/12.0F, 1E-4F );    
    dscan_convexpolyhedron_destroy( tetra );
}

void Check_convex( bool coverage )
{
    if( coverage )
    {
        Check_convex_check();
        Check_convex_3D_check();
    }
    Check_convex_intersect();
    Check_convex_project();
    Check_convex_3D();
}

void CheckPath
(
 float* data,
 int sides,
 vector2_t start,
 vector2_t end
)
{
    vector<vector2_t> polyPoints;
    for( int i = 0; i < sides; ++i )
    {
        polyPoints.push_back( MakeVector( data[2*i], data[2*i + 1] ) );
    }
    convexpolygon_t poly = MakeConvexPolygon( polyPoints );

    dscan_lightpath_t lp = dscan_lightpath_create( start, end );
    dscan_material_t external = dscan_material_create( "air", 3 );
    dscan_material_t internal = dscan_material_create( "crownimpure", 11 );
    int rc;
    while( (rc = dscan_convexpolygon_trace( lp, external, internal, &poly, 0.0001F, 1E-4F ) ) );

    int p = dscan_lightpath_points( lp );
    int e = dscan_lightpath_elements( lp );

    dscan_vector2_t* points = new dscan_vector2_t[p];
    int* pointCnts = new int[e];
    int* types = new int[e];
    float* weights = new float[e];
    dscan_lightpath_marshall( lp, 0, e, points, pointCnts, types, weights );
    
    int pi = 0;
    for( int ei = 0; ei < e; ++ei )
    {
        cout << "Type " << types[ei] << " "
             << "Weight " << weights[ei] << " "
             << "Points";
        for( int i = 0; i < pointCnts[ei]; ++i, ++pi )
        {
            cout << " <" << points[pi].x << "," << points[pi].y << ">";
        }
        cout << endl;
    }

    delete[] weights;
    delete[] types;
    delete[] pointCnts;
    delete[] points;
    dscan_material_destroy( external );
    dscan_material_destroy( internal );
    dscan_lightpath_destroy( lp );
    convexpolygon_destroy( poly );
}

void Check_path()
{
    float p1[] = {-1,0, 0,-1, 1,0, 0,1};
    vector2_t s1 = MakeVector( 0.5F, 2 );
    vector2_t e1 = MakeVector( -0.5F, 2 );
    CheckPath( p1, 4, s1, e1 );

    dscan_roundcut_measurements_t m;
    dscan_roundcut_errors_t err;
    m.table_pct = 0.57F;
    m.crownheight_pct = 0.15F; /*extra*/
    m.pavilionheight_pct = 0.445F; /*extra*/
    m.girdle_pct = 0.035F; /*extra*/
    m.star_pct = 0.55F;
    m.pavilion_pct = 0.75F;
    m.culet_pct = 0.0F;
    m.girdleedges = 4;
    dscan_convexpolyhedron_t poly = dscan_roundcut_create( m, &err );
    dscan_material_t external = dscan_material_create( "air", 3 );
    dscan_material_t internal = dscan_material_create( "diamond", 7 );
    dscan_convexpolygon_t cross = dscan_roundcut_crosssection( poly, 0 );
    int mx = 0;
    while(mx < 4500)
    {
        matrix3_t rot = RandomTransform();
        vector2_t start = NonHMultiply( rot, MakeVector( 2, 2 ) );
        vector2_t end = NonHMultiply( rot, MakeVector( -2, 2 ) );
        dscan_lightpath_t lp = dscan_lightpath_create( start, end );
        while( dscan_convexpolygon_trace( lp, external, internal, cross, 1E-4F, 1E-4F ) );
        int elementsCnt = dscan_lightpath_elements( lp );
        if( elementsCnt > mx )
        {
            mx = elementsCnt;
            cout << mx << endl;
        }
        dscan_lightpath_destroy( lp );
    }
    dscan_convexpolyhedron_destroy( poly );
    dscan_material_destroy( external );
    dscan_material_destroy( internal );
    dscan_convexpolygon_destroy( cross );
}

int main( int argc, char** argv )
{
//    bool coverage = 0;
    if( argc == 2 && string( argv[1] ) == "COVERAGE" )
    {
        CVXITERS = 1;
//        coverage = 1;
    }

    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW );

//    Check_vec();
//    Check_cut();
//    Check_convex( coverage );
    Check_path();

    cout << "SUCCESS" << endl;
}
