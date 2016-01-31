#define _GNU_SOURCE
#include <fenv.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dscan.h"

static const int ITERS = 10;

float* make_float( float x )
{
    float* retval = malloc( sizeof(float) );
    *retval = x;
    return retval;
}

int* make_int( int x )
{
    int* retval = malloc( sizeof(int) );
    *retval = x;
    return retval;
}

float to_radians( float d )
{
    return d * 3.14159265F / 180.0F;
}

int main( int argc, char** argv )
{
    dscan_roundcut_input_t in,out;
    dscan_roundcut_measurements_t m;
    dscan_roundcut_errors_t err;
    dscan_convexpolyhedron_t poly;
    dscan_vector2_t source_start, source_end;
    dscan_convexpolygon_t cross_section;
    dscan_material_t internal, external;
    dscan_observer_t observer;
    dscan_lightpath_t path;
    dscan_matrix3_t display_transform;
    int pointsCnt;
    int elementsCnt;
    dscan_vector2_t* points;
    int* pointCnts;
    int* types;
    float* weights;
    int* buffer = malloc( sizeof(int) * 480 * 618 );
    int i, j, iter;
    time_t start, end;
    /*float theta;*/

    feenableexcept( FE_DIVBYZERO | /*FE_INVALID |*/ FE_OVERFLOW |
                    FE_UNDERFLOW );

    memset( &in, 0, sizeof(dscan_roundcut_input_t) );
    memset( &out, 0, sizeof(dscan_roundcut_input_t) );
    in.width_min = make_float(0.0690F);
    in.width_max = make_float(0.0694F);
    in.height = make_float(0.0426F);
    in.table_pct = make_float(0.57F);
    in.crown_angle = make_float(to_radians(34.5F));
    in.pavilion_angle = make_float(to_radians(40.8F));
    in.culet_pct = make_float(0.0F);
    in.star_pct = make_float(0.55F);
    in.pavilion_pct = make_float(0.75F);
    in.girdleedges = make_int(4);
    dscan_roundcut_measure( in, &out );

    /*m.min_width = 0.0690F;*/
    /*m.max_width = 0.0694F;*/
    /*m.height = 0.0426F;*/
    /*m.table_pct = 0.57F;*/
    /*m.height_pct = 0.615F;*/ /*extra*/
    /*m.crownheight_pct = 0.15F;*/ /*extra*/
    /*m.pavilionheight_pct = 0.445F;*/ /*extra*/
    /*m.crownangle_deg = 34.5F;*/
    /*m.pavilionangle_deg = 40.8F;*/
    /*m.girdle_pct = 0.035F;*/ /*extra*/
    /*m.star_pct = 0.55F;*/
    /*m.pavilion_pct = 0.75F;*/
    /*m.culet_pct = 0.0F;*/
    /*m.girdleedges = 4;*/
    m = dscan_roundcut_complete( in, out );
    poly = dscan_roundcut_create( m, &err );
    /* dscan_convexpolyhedron_check( poly ); */

    internal = dscan_material_create( "crownimpure", 11 );
    external = dscan_material_create( "air", 3 );
    observer = dscan_observer_create( "default", 7 );
    cross_section = dscan_roundcut_crosssection( poly, 0 );
    source_start.x = 0.0F;
    source_start.y = 0.5F;
    /*source_end.x = -(0.57F / 2.0F);
      source_end.y = 0.5F;*/
    /*source_end.x = -0.5F;
      source_end.y = 0.5F;*/
    /*source_end.x = -0.36249995F;
      source_end.y =  0.65937496F;*/
    source_end.x = -0.67499995F;
    source_end.y =  0.42812496F;
    display_transform.d[0] = 320.0F;
    display_transform.d[1] = 0.0F;
    display_transform.d[2] = 240.0F;
    display_transform.d[3] = 0.0F;
    display_transform.d[4] = 320.0F;
    display_transform.d[5] = 309.0F;
    display_transform.d[6] = 0.0F;
    display_transform.d[7] = 0.0F;
    display_transform.d[8] = 1.0F;
    /*for( int theta = -0.25F; theta < 0.25F; theta += 0.01F )*/
    start = time(0);
    for( iter = 0; iter < ITERS; ++iter )
    {
/*        source_end.x = source_start.x - 0.5F * cos(theta);
          source_end.y = source_start.y + 0.5F * sin(theta);*/

        path = dscan_lightpath_create(source_start, source_end);
        for( i = 0; i < 10; ++i )
        {
            int cnt = dscan_convexpolygon_trace( path, external, internal,
                                                 cross_section, 1E-4F, 1E-4F );
            if( cnt <= 0 ) break;
        }
        pointsCnt = dscan_lightpath_points( path );
        elementsCnt = dscan_lightpath_elements( path );
        points = malloc( sizeof(dscan_vector2_t) * pointsCnt );
        pointCnts = malloc( sizeof(int) * elementsCnt );
        types = malloc( sizeof(int) * elementsCnt );
        weights = malloc( sizeof(float) * elementsCnt );
        dscan_lightpath_marshall( path, 0, elementsCnt, points, pointCnts, 
                                  types, weights );
        memset( buffer, 0, 480*618 * sizeof(int) );
        for( i = 0, j = 0; i < elementsCnt; ++i )
        {
            if( 1 == types[i] )
            {
                dscan_vector2_t start = points[j];
                dscan_vector2_t end = points[j+1];
                dscan_vector2_t ray = points[j+2];
                float weight = weights[i];
                dscan_prism_t prism = dscan_prism_create( start, end, ray,
                                                          internal, external,
                                                          observer, 
                                                          display_transform,
                                                          weight );
                dscan_prism_rasterize( prism, buffer, 0, 480, 480, 618 );
                dscan_prism_destroy( prism );
            }
            j += pointCnts[i];
        }
        dscan_lightpath_destroy( path );

        free( points );
        free( pointCnts );
        free( types );
        free( weights );
    }
    end = time(0);
    fprintf( stderr, "%dms %d elements\n", (int)((end-start)*1000/ITERS), elementsCnt );

    dscan_material_destroy( internal );
    dscan_material_destroy( external );
    dscan_observer_destroy( observer );
    dscan_convexpolygon_destroy( cross_section );
    dscan_convexpolyhedron_destroy( poly );

    printf( "480 618\n" );
    for( i = 0; i < 480*618; ++i ) {
        int argb = buffer[i];
        printf( "%d %d %d %d\n",
                (argb >> 24) & 0xFF,
                (argb >> 16) & 0xFF,
                (argb >> 8) & 0xFF,
                argb & 0xFF );
    }

    free( buffer );

    return 0;
}
