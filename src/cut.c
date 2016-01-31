#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "dscan.h"
#include "dscanimp.h"

#define BT 0
#define BR 1
#define BB 2
#define U0R 3
#define PT 4
#define PR 5
#define L0R 6
/*For use with culet*/
#define PB 7 

typedef struct plane_t
{
    vector3_t ref;
    vector3_t norm;
} plane_t;

static vector3_t make_vector3( float x, float y, float z )
{
    vector3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static vector2_t make_vector2(float x_, float y_) 
{
    vector2_t v;
    v.x = x_;
    v.y = y_;
    return v;
}

static edge_t make_edge(int a, int b)
{
    edge_t retval;
    retval.a = a;
    retval.b = b;
    return retval;
}

static face_t alloc_face3(int a, bool ra, int b, bool rb, int c)
{
    face_t retval;
    retval.edges = malloc(sizeof(int) * 3);
    retval.edgesCnt = 3;
    retval.edges[0] = a;
    retval.edges[1] = b;
    retval.edges[2] = c;

    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static face_t alloc_face4(int a, bool ra, int b, bool rb, int c, int d)
{
    face_t retval;
    retval.edges = malloc(sizeof(int) * 4);
    retval.edgesCnt = 4;
    retval.edges[0] = a;
    retval.edges[1] = b;
    retval.edges[2] = c;
    retval.edges[3] = d;
    
    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static face_t alloc_face5(int a, bool ra, int b, bool rb, int c, int d, int e)
{
    face_t retval;
    retval.edges = malloc(sizeof(int) * 5);
    retval.edgesCnt = 5;
    retval.edges[0] = a;
    retval.edges[1] = b;
    retval.edges[2] = c;
    retval.edges[3] = d;
    retval.edges[4] = e;
    
    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static face_t alloc_face6(int a, bool ra, int b, bool rb, int c, int d, int e, int f)
{
    face_t retval;
    retval.edges = malloc(sizeof(int) * 6);
    retval.edgesCnt = 6;
    retval.edges[0] = a;
    retval.edges[1] = b;
    retval.edges[2] = c;
    retval.edges[3] = d;
    retval.edges[4] = e;
    retval.edges[5] = f;

    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static face_t alloc_face_n1(int a, bool ra, int b, bool rb, int* r, int len) 
{
    face_t retval;
    int i;
    retval.edges = malloc(sizeof(int) * (2 + len));
    retval.edgesCnt = 2 + len;
    retval.edges[0] = a;
    retval.edges[1] = b;
    for(i = 0; i < len; ++i) 
    {
        retval.edges[i + 2] = r[i];
    }

    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static face_t make_face_n2(int* r, int len, bool ra, bool rb) 
{
    face_t retval;
    retval.edges = r;
    retval.edgesCnt = len;
    retval.reverse0 = ra;
    retval.reverse1 = rb;
    return retval;
}

static matrix3_t make_yrotation(int i)
{
    matrix3_t r;
    double h = sqrt(2.0) / 2.0;
    r.d[1] = r.d[3] = r.d[5] = r.d[7] = 0.0F;
    r.d[4] = 1.0F;
    switch(i) {
    default:
        assert(0);
    case 0:
        r.d[0] = r.d[8] = 1.0F;
        r.d[2] = r.d[6] = 0.0F;
        break;
    case 1:
        r.d[0] = r.d[2] = r.d[8] = h;
        r.d[6] = -h;
        break;
    case 2:
        r.d[0] = r.d[8] = 0.0F;
        r.d[2] = 1.0F;
        r.d[6] = -1.0F;
        break;
    case 3:
        r.d[2] = h;
        r.d[0] = r.d[6] = r.d[8] = -h;
        break;
    case 4:
        r.d[0] = r.d[8] = -1.0F;
        r.d[2] = r.d[6] = 0.0F;
        break;
    case 5:
        r.d[6] = h;
        r.d[0] = r.d[2] = r.d[8] = -h;
        break;
    case 6:
        r.d[0] = r.d[8] = 0.0F;
        r.d[2] = -1.0F;
        r.d[6] = 1.0F;
        break;
    case 7:
        r.d[0] = r.d[6] = r.d[8] = h;
        r.d[2] = -h;
        break;
    }
    
    return r;
}

static vector3_t intersect( plane_t* plane, vector3_t iref, vector3_t idir ) 
{
    vector3_t rel = minus3( iref, plane->ref );
    vector3_t i = intersect_plane( plane->norm, rel, idir );
    i = add3( i, plane->ref );
    return i;
}

static plane_t make_plane
(
 vector3_t p1, 
 vector3_t p2, 
 vector3_t p3
)
{
    plane_t retval;
    vector3_t l0 = minus3( p2, p1 );
    vector3_t l1 = minus3( p3, p1 );
    retval.ref = p1;
    retval.norm = cross3(l0, l1);
    retval.norm = normalize3( retval.norm );
    return retval;
}

static void copy_down
(
 vector2_t* s,
 vector3_t t
)
{
    s->x = t.x;
    s->y = t.y;
}

static vector3_t mid3( vector3_t a, vector3_t b )
{
    vector3_t retval;
    retval.x = a.x / 2.0F + b.x / 2.0F;
    retval.y = a.y / 2.0F + b.y / 2.0F;
    retval.z = a.z / 2.0F + b.z / 2.0F;
    return retval;
}

#define TABLE_PCT 0
#define CROWN_PCT 1
#define PAVILION_PCT 2
#define GIRDLE_PCT 3
#define CULET_PCT 4
#define WIDTH_MIN 5
#define WIDTH_MAX 6
#define HEIGHT 7
#define PAVILION_ANGLE 8
#define WIDTH_AVG 9
#define HEIGHT_PCT 10
#define CROWN_ANGLE 11

static float** raw_field_input(dscan_roundcut_input_t* in, int field)
{
    switch(field)
    {
    case TABLE_PCT: return &in->table_pct;
    case CROWN_PCT: return &in->crownheight_pct;
    case PAVILION_PCT: return &in->pavilionheight_pct;
    case GIRDLE_PCT: return &in->girdle_pct;
    case CULET_PCT: return &in->culet_pct;
    case WIDTH_MIN: return &in->width_min;
    case WIDTH_MAX: return &in->width_max;
    case HEIGHT: return &in->height;
    case PAVILION_ANGLE: return &in->pavilion_angle;
    case WIDTH_AVG: return &in->width_avg;
    case HEIGHT_PCT: return &in->height_pct;
    case CROWN_ANGLE: return &in->crown_angle;
    default: assert(0); return 0;
    }   
}

static float* get_field_input(dscan_roundcut_input_t* in, int field)
{
    return *raw_field_input(in, field);
}

static int has_field_input(dscan_roundcut_input_t* in, int field)
{
    return !!get_field_input(in, field);
}

static float* make_field_input(dscan_roundcut_input_t* in, int field)
{
    float** ptr = raw_field_input(in, field);
    *ptr = malloc(sizeof(float));
    return *ptr;
}

static float* get_field_meas(dscan_roundcut_measurements_t* in, int field)
{
    switch(field)
    {
    case TABLE_PCT: return &in->table_pct;
    case CROWN_PCT: return &in->crownheight_pct;
    case PAVILION_PCT: return &in->pavilionheight_pct;
    case GIRDLE_PCT: return &in->girdle_pct;
    case CULET_PCT: return &in->culet_pct;
    default: assert(0); return 0;
    }   
}

/*
static int has_field_meas(dscan_roundcut_measurements_t* in, int field)
{
    return !!get_field_meas(in, field);
}
*/

static int should_compute
(
 dscan_roundcut_input_t* in,
 dscan_roundcut_input_t* out,
 const int fields[],
 float data[],
 float** result,
 int n,
 int priority
)
{
    int i;

    if( (!priority && has_field_input( in, fields[n-1] )) ||
        has_field_input( out, fields[n-1] ) )
    {
        return 0;
    }

    for( i = 0; i < n-1; ++i )
    {
        if( has_field_input( out, fields[i] ) )
        {
            data[i] = *get_field_input( out, fields[i] );
        }
        else if( has_field_input( in, fields[i] ) )
        {
            data[i] = *get_field_input( in, fields[i] );
        }
        else
        {
            return -1;
        }
    }

    *result = get_field_input( out, fields[n-1] );

    if( !*result )
    {
        *result = make_field_input( out, fields[n-1] );
    }

    return 1;
}

#define FUNC2(NM,F0,F1,R,C,P)                                           \
int NM(dscan_roundcut_input_t* in, dscan_roundcut_input_t* out)         \
{                                                                       \
const int Fields[] = {F0,F1,R};                                         \
    float d[2], *r;                                                     \
    int retval;                                                         \
    if( 1 == (retval = should_compute(in, out, Fields, d, &r, 3, P )) ) \
    {                                                                   \
        *r = C;                                                         \
    }                                                                   \
    return retval;                                                      \
}

#define FUNC3(NM,F0,F1,F2,R,C,P)                                        \
int NM(dscan_roundcut_input_t* in, dscan_roundcut_input_t* out)         \
{                                                                       \
    const int Fields[] = {F0,F1,F2,R};                                  \
    float d[3], *r;                                                     \
    int retval;                                                         \
    if( 1 == (retval = should_compute(in, out, Fields, d, &r, 4, P )) ) \
    {                                                                   \
        *r = C;                                                         \
    }                                                                   \
    return retval;                                                      \
}

FUNC2(rcm_www_w,WIDTH_MIN,WIDTH_MAX,WIDTH_AVG,d[0]/2.0f+d[1]/2.0f,1)
FUNC2(rcm_whp_w,HEIGHT,HEIGHT_PCT,WIDTH_AVG,d[0]/d[1],0)
FUNC2(rcm_whp_h,HEIGHT_PCT,WIDTH_AVG,HEIGHT,d[0]*d[1],0)
FUNC2(rcm_whp_p,HEIGHT,WIDTH_AVG,HEIGHT_PCT,d[0]/d[1],1)
FUNC2(rcm_tcp_t,CROWN_ANGLE,CROWN_PCT,TABLE_PCT,1-d[1]/(2*tan(d[0])),0)
FUNC2(rcm_tcp_c,CROWN_PCT,TABLE_PCT,CROWN_ANGLE,atan(2*d[0]/(1-d[1])),0)
FUNC2(rcm_tcp_p,TABLE_PCT,CROWN_ANGLE,CROWN_PCT,tan(d[1])*(1-d[0])/2.0F,1)
FUNC3(rcm_hcpg_h,CROWN_PCT,PAVILION_PCT,GIRDLE_PCT,HEIGHT_PCT,d[0]+d[1]+d[2],0)
FUNC3(rcm_hcpg_c,PAVILION_PCT,GIRDLE_PCT,HEIGHT_PCT,CROWN_PCT,d[2]-d[0]-d[1],0)
FUNC3(rcm_hcpg_p,GIRDLE_PCT,HEIGHT_PCT,CROWN_PCT,PAVILION_PCT,d[1]-d[0]-d[2],0)
FUNC3(rcm_hcpg_g,HEIGHT_PCT,CROWN_PCT,PAVILION_PCT,GIRDLE_PCT,d[0]-d[1]-d[2],1)
FUNC2(rcm_pac_a,CULET_PCT,PAVILION_PCT,PAVILION_ANGLE,atan(2*d[1]/(1-d[0])),0)
FUNC2(rcm_pac_c,PAVILION_PCT,PAVILION_ANGLE,CULET_PCT,1-d[0]/(2*tan(d[1])),0)
FUNC2(rcm_pac_p,PAVILION_ANGLE,CULET_PCT,PAVILION_PCT,tan(d[0])*(1-d[1])/2.0F,1)

int dscan_roundcut_measure(dscan_roundcut_input_t in, dscan_roundcut_input_t* out)
{
    int (*funcs[14])(dscan_roundcut_input_t*,dscan_roundcut_input_t*);
    int tail = 13;
    int advanced = 1;
    funcs[0] = &rcm_www_w;
    funcs[1] = &rcm_whp_w;
    funcs[2] = &rcm_whp_h;
    funcs[3] = &rcm_whp_p;
    funcs[4] = &rcm_tcp_t;
    funcs[5] = &rcm_tcp_c;
    funcs[6] = &rcm_tcp_p;
    funcs[7] = &rcm_hcpg_h;
    funcs[8] = &rcm_hcpg_c;
    funcs[9] = &rcm_hcpg_p;
    funcs[10] = &rcm_hcpg_g;
    funcs[11] = &rcm_pac_p;
    funcs[12] = &rcm_pac_a;
    funcs[13] = &rcm_pac_c;

    while( advanced && tail >= 0 )
    {
        int i = 0;
        advanced = 0;
        while( i <= tail )
        {
            int rc = funcs[i](&in, out);
            if( 1 == rc )
            {
                advanced = 1;
                ++i;
            }
            else if( -1 == rc )
            {
                ++i;
            }
            else
            {
                funcs[i] = funcs[tail--];
            }
        }
    }

    return tail < 0 && in.star_pct && in.pavilion_pct && in.girdleedges;
}
 
dscan_roundcut_measurements_t dscan_roundcut_complete(dscan_roundcut_input_t given, dscan_roundcut_input_t calculated)
{
    dscan_roundcut_measurements_t retval;
    int i;
    for( i = 0; i < 5; ++i )
    {
        if( has_field_input( &calculated, i ) )
        {
            *get_field_meas( &retval, i ) = *get_field_input( &calculated, i );
        }
        else
        {
            assert( has_field_input( &given, i ) );
            *get_field_meas( &retval, i ) = *get_field_input( &given, i );
        }
    }
    if( given.star_pct )
    {
        retval.star_pct = *given.star_pct;
    }
    if( given.pavilion_pct )
    {
        retval.pavilion_pct = *given.pavilion_pct;
    }
    if( given.girdleedges )
    {
        retval.girdleedges = *given.girdleedges;
    }
    return retval;
}

convexpolyhedron_t* dscan_roundcut_create
(
 roundcut_measurements_t m,
 roundcut_errors_t* err
)
{
    /*
     * Naming Convention:
     * 
     * {facet}{position}
     * 
     * Facets: s : Star facets - Triangle connected to top table. b : Bezel
     * facets - Kite shapes, upper formed by star facets. u# : Upper girdle
     * - Triangles connected to upper girdle. l# : Lower girdle - Triangles
     * connected to lower girdle. p : Pavilion - Kite shapes, connected to
     * culet point.
     * 
     * Positions: t : Top point - For use with bezel, upper girdle and
     * pavilion. l : Left point - For use with all. i# : Internal points -
     * For use with girdle connections. r : Right point - For use with all.
     * b : Bottom point - For use with bezel, lower facets.
     */

    convexpolyhedron_t* retval;
    int i, j, k;
    
    float tableh = m.table_pct / 2.0F;
    float crown = m.crownheight_pct;
    float star = m.star_pct;
    float pavilion = m.pavilion_pct;
    float girdle = m.girdle_pct;
    int girdleedges = m.girdleedges;
    float depth = girdle + m.pavilionheight_pct;
    float culet = m.culet_pct;
    float thetad = (float) (2.0 * M_PI / 8.0);
    bool has_culet = (m.culet_pct > 0.0F);

    const int ipoints = 4 + 3 + 4 * girdleedges + (has_culet ? 1 : 0);
    const  int npoints = 8 * ipoints + (has_culet ? 0 : 1);
    vector3_t* points = malloc(sizeof(vector3_t)*npoints);
    const int iedges = 6 + 4 + 4 * (girdleedges + 1) + 2 * girdleedges + 
        (has_culet ? 1 : 0);
    const int nedges = 8 * iedges;
    edge_t* edges = malloc(sizeof(edge_t)*nedges);
    const int ifaces = 4 + 3 + 2 * girdleedges;
    const int nfaces = 8 * ifaces + 1 + (has_culet ? 1 : 0);
    face_t* faces = malloc(sizeof(face_t)*nfaces);
    /* E = 8*(6+4+4(g+1)+2g+c) = 112+48g+8c
     * F = 8*(4+3+2g)+1+c = 57+16g+c
     * V = 8*(4+3+4g+c)+(1-c) = 57+32g+7c
     * F+V-2 = 112+48g+8c
     */

    vector3_t* lp = malloc(sizeof(vector3_t)*ipoints);

    const int U0Iref = (has_culet ? PB+1 : L0R+1 );
    const int U1Iref = U0Iref + girdleedges;
    const int L0Iref = U1Iref + girdleedges;
    const int L1Iref = L0Iref + girdleedges;
    const int CL = (has_culet ? ipoints - 1 : -1);

    /* Compute primary vectors */
    vector2_t a0 = make_vector2(1.0F, 0.0F);
    vector2_t a1 = make_vector2(cos(thetad / 2.0), sin(thetad / 2.0));
    vector2_t a2 = make_vector2(cos(thetad), sin(thetad));

    float girdleRadius = 0.5F / cos(thetad / (girdleedges * 4));
    bool has_error = false;

    /* Point work */
    vector3_t u1r;
    vector3_t l1r;
    vector3_t pb;

    vector3_t yAxis;

    plane_t b;
    float sScale;
    plane_t p;
    float pScale;
    plane_t u0, u1, l0, l1;

    /* edges */
    const int BT_BTp = 0;
    const int BT_BL = 1;
    const int BL_BB = 2;
    const int BB_BR = 3;
    const int BR_BT = 4;
    const int BR_U0R = 5;
    const int PT_PL = 6;
    const int PB_PR = 7;
    const int PR_PT = 8;
    const int PR_L0R = 9;
    const int U0Ieref = 10;
    const int U1Ieref = U0Ieref + girdleedges + 1;
    const int L0Ieref = U1Ieref + girdleedges + 1;
    const int L1Ieref = L0Ieref + girdleedges + 1;
    const int U0I_L0Iref = L1Ieref + girdleedges + 1;
    const int U1I_L1Iref = U0I_L0Iref + girdleedges;
    const int PB_C = U1I_L1Iref + girdleedges;

    edge_t* le = malloc(sizeof(edge_t) * iedges);

    /* faces */
    const int B = 0;
    const int S = 1;
    const int U0 = 2;
    const int U1 = 3;
    const int P = 4;
    const int L0 = 5;
    const int L1 = 6;
    const int G0ref = 7;
    const int G1ref = G0ref + girdleedges;

    face_t* lf = malloc(sizeof(face_t) * ifaces);

    int* fu0;
    int* fu1;
    int* fl0;
    int* fl1;
    int* table;

    /* initialize errs */
    err->min_girdle_pct = girdle;

    /* Compute major markers. */
    lp[BT] = make_vector3( tableh * a0.x, crown, tableh * a0.y );
    lp[BB] = make_vector3( 0.5F * a0.x, 0.0F, 0.5F * a0.y );
    lp[U0R] = make_vector3( 0.5F * a1.x, 0.0F, -0.5F * a1.y );
    u1r = make_vector3( 0.5F * a2.x, 0.0F, -0.5F * a2.y );
    lp[PT] = make_vector3( 0.5F * a0.x, -girdle, 0.5F * a0.y );
    lp[L0R] = make_vector3( 0.5F * a1.x, -girdle, -0.5F * a1.y );
    l1r = make_vector3( 0.5F * a2.x, -girdle, -0.5F * a2.y );
    pb = make_vector3( 0.0F, -depth, 0.0F );

    /* Compute intersect vector */
    yAxis = make_vector3( 0.0F, 1.0F, 0.0F );

    /* Compute bezel facet sides. */
    b = make_plane(lp[BT], lp[BB], make_vector3( lp[BT].x, lp[BT].y, -1.0F ) );
    sScale = tableh * a1.x + (0.5F - tableh * a1.x) * star;
    lp[BR] = intersect( &b, make_vector3( sScale * a1.x, 0, -sScale * a1.y ),
                        yAxis );

    /* Compute pavilion facet sides. */
    p = make_plane( pb, lp[PT], make_vector3( pb.x, pb.y, -1.0F ) );
    pScale = 0.5F * (1.0F - pavilion);
    lp[PR] = intersect( &p, make_vector3( pScale * a1.x, 0, -pScale * a1.y ),
                        yAxis );

    /* Compute culet. */
    if( has_culet )
    {
        float culet_depth = girdle + (depth-girdle) * (1-culet/0.5F);
        plane_t c;
        vector3_t l = minus3( pb, lp[PR] );
        c.ref = make_vector3( 0, -culet_depth, 0 );
        c.norm = make_vector3( 0, -1.0F, 0 );
        lp[PB] = intersect( &c, lp[PR], l );
    }

    /* Compute girdle. */
    u0 = make_plane(lp[BR], lp[BB], lp[U0R]);
    u1 = make_plane(lp[BR], lp[U0R], u1r);
    l0 = make_plane(lp[PR], lp[PT], lp[L0R]);
    l1 = make_plane(lp[PR], lp[L0R], l1r);
    for (i = 0; i < girdleedges; ++i) {
        float div = 4 * girdleedges; /*16.0F; TODO*/
        int incr = 2 * girdleedges; /*2; TODO*/
        vector3_t t0 = make_vector3( 
         girdleRadius * cos(thetad / div * (2 * i + 1) ), 
         0.0F, -girdleRadius * sin(thetad / div * (2 * i + 1)) );
        vector3_t t1 = make_vector3(
         girdleRadius * cos(thetad / div * (incr + 2 * i + 1)), 
         0.0F, -girdleRadius * sin(thetad / div * (incr + 2 * i + 1)) );
        lp[U0Iref + i] = intersect(&u0, t0, yAxis);
        lp[U1Iref + i] = intersect(&u1, t1, yAxis);
        lp[L0Iref + i] = intersect(&l0, t0, yAxis);
        lp[L1Iref + i] = intersect(&l1, t1, yAxis);

        /* Sync */
        lp[L0Iref + i].x = lp[U0Iref + i].x;
        lp[L0Iref + i].z = lp[U0Iref + i].z;
        lp[L1Iref + i].x = lp[U1Iref + i].x;
        lp[L1Iref + i].z = lp[U1Iref + i].z;

        /* Check */
        if( lp[U0Iref + i].y <= lp[L0Iref + i].y )
        {
            float mn = girdle + (lp[L0Iref + i].y - lp[U0Iref + i].y);
            if( mn > err->min_girdle_pct )
            {
                err->min_girdle_pct = mn;
            }
            has_error = true;
        }
        if( lp[U1Iref + i].y <= lp[L1Iref + i].y )
        {
            float mn = girdle + (lp[L1Iref + i].y - lp[U1Iref + i].y);
            if( mn > err->min_girdle_pct )
            {
                err->min_girdle_pct = mn;
            }
            has_error = true;
        }
    }

    /* Dupliate points */
    for( i = 0; i < 8; ++i )
    {
        matrix3_t r = make_yrotation(i);
        for( j = 0; j < ipoints; ++j )
        {
            vector3_t* p = points + i*ipoints + j;
            vector3_t* src = lp + j;
            float x = src->x;
            float y = src->y;
            float z = src->z;
            p->x = r.d[0] * x + r.d[1] * y + r.d[2] * z;
            p->y = r.d[3] * x + r.d[4] * y + r.d[5] * z;
            p->z = r.d[6] * x + r.d[7] * y + r.d[8] * z;
        }
    }
    if( !has_culet )
    {
        points[8 * ipoints] = pb;
    }
    free(lp);

    /* Place edges */

    le[BT_BTp] = make_edge(BT, BT + ipoints);
    le[BT_BL] = make_edge(BT, BR - ipoints);
    le[BL_BB] = make_edge(BR - ipoints, BB);
    le[BB_BR] = make_edge(BB, BR);
    le[BR_BT] = make_edge(BR, BT);
    le[BR_U0R] = make_edge(BR, U0R);
    le[PT_PL] = make_edge(PT, PR - ipoints);
    le[PB_PR] = make_edge((has_culet ? PB : CL), PR);
    le[PR_PT] = make_edge(PR, PT);
    le[PR_L0R] = make_edge(PR, L0R);
    le[U0Ieref + 0] = make_edge(BB, U0Iref + 0);
    le[U1Ieref + 0] = make_edge(U0R, U1Iref + 0);
    le[L0Ieref + 0] = make_edge(PT, L0Iref + 0);
    le[L1Ieref + 0] = make_edge(L0R, L1Iref + 0);
    for( i = 0; i < girdleedges - 1; ++i )
    {
        le[U0Ieref + 1 + i] = make_edge(U0Iref + i, U0Iref + i + 1);
        le[U1Ieref + 1 + i] = make_edge(U1Iref + i, U1Iref + i + 1);
        le[L0Ieref + 1 + i] = make_edge(L0Iref + i, L0Iref + i + 1);
        le[L1Ieref + 1 + i] = make_edge(L1Iref + i, L1Iref + i + 1);
    }
    le[U0Ieref + girdleedges] = make_edge(U0Iref+girdleedges-1, U0R);
    le[U1Ieref + girdleedges] = make_edge(U1Iref+girdleedges-1, BB + ipoints);
    le[L0Ieref + girdleedges] = make_edge(L0Iref+girdleedges-1, L0R);
    le[L1Ieref + girdleedges] = make_edge(L1Iref+girdleedges-1, PT + ipoints);
    for( i = 0; i < girdleedges; ++i )
    {
        le[U0I_L0Iref + i] = make_edge(U0Iref + i, L0Iref + i);
        le[U1I_L1Iref + i] = make_edge(U1Iref + i, L1Iref + i);
    }
    if( has_culet )
    {
        le[PB_C] = make_edge( PB, PB - ipoints );
    }

    /* Duplicate edges */
    for( i = 0; i < 8; ++i )
    {
        for( j = 0; j < iedges; ++j )
        {
            edge_t* edge = edges + i*iedges + j;
            edge->a = (le[j].a + (8+i)*ipoints) % (8*ipoints); 
            edge->b = (le[j].b + (8+i)*ipoints) % (8*ipoints);
        }
        if( !has_culet )
        {
            edges[i * iedges + PB_PR].a = npoints - 1;
        }
    }
    free(le);

    /* Place faces */
    lf[B] = alloc_face4(BT_BL, 0, BL_BB, 0, BB_BR, BR_BT);
    lf[S] = alloc_face3(BR_BT, 1, BT_BL + iedges, 1, BT_BTp);
    if( has_culet )
    {
        lf[P] = alloc_face5(PR_PT, 0, PT_PL, 0, PB_PR - iedges, PB_C, PB_PR );
    }
    else
    {
        lf[P] = alloc_face4(PR_PT, 0, PT_PL, 0, PB_PR - iedges, PB_PR);
    }
    fu0 = malloc(sizeof(int) * (girdleedges + 1));
    fu1 = malloc(sizeof(int) * (girdleedges + 1));
    fl0 = malloc(sizeof(int) * (girdleedges + 1));
    fl1 = malloc(sizeof(int) * (girdleedges + 1));
    for (i = 0; i < girdleedges + 1; ++i) {
        fu0[i] = U0Ieref + i;
        fu1[i] = U1Ieref + i;
        fl0[girdleedges - i] = L0Ieref + i;
        fl1[girdleedges - i] = L1Ieref + i;
    }
    lf[U0] = alloc_face_n1(BR_U0R, 1, BB_BR, 1, fu0, girdleedges + 1);
    lf[U1] = alloc_face_n1(BL_BB + iedges, 1, BR_U0R, 0, fu1, girdleedges + 1);
    lf[L0] = alloc_face_n1(PR_PT, 1, PR_L0R, 0, fl0, girdleedges + 1);
    lf[L1] = alloc_face_n1(PR_L0R, 1, PT_PL + iedges, 1, fl1, girdleedges + 1);
    lf[G0ref + 0] = alloc_face6(L0Ieref + 0, 0, U0I_L0Iref + 0, 1, 
                                U0Ieref + 0, U1Ieref + girdleedges - iedges, 
                                U1I_L1Iref + girdleedges - 1 - iedges, 
                                L1Ieref + girdleedges - iedges);
    lf[G1ref + 0] = alloc_face6(L1Ieref + 0, 0, U1I_L1Iref + 0, 1, 
                                U1Ieref + 0, U0Ieref + girdleedges, 
                                U0I_L0Iref + girdleedges - 1, 
                                L0Ieref + girdleedges);
    for( i = 1; i < girdleedges; ++i )
    {
        lf[G0ref + i] = alloc_face4(L0Ieref + i, 0, U0I_L0Iref + i, 1, 
                                    U0Ieref + i, U0I_L0Iref + i - 1);
        lf[G1ref + i] = alloc_face4(L1Ieref + i, 0, U1I_L1Iref + i, 1, 
                                    U1Ieref + i, U1I_L1Iref + i - 1);
    }
    free(fu0);
    free(fu1);
    free(fl0);
    free(fl1);

    /* Duplicate faces */
    for( i = 0; i < 8; ++i )
    {
        for( j = 0; j < ifaces; ++j )
        {
            face_t* face;
            int len = lf[j].edgesCnt;
            int* r = malloc(sizeof(int) * len);
            for( k = 0; k < len; ++k )
            {
                r[k] = (lf[j].edges[k] + (8 + i) * iedges) % (8 * iedges);
            }
            face = faces + i*ifaces + j;
            face->edges = r;
            face->edgesCnt = len;
            face->reverse0 = lf[j].reverse0;
            face->reverse1 = lf[j].reverse1;
        }
    }

    /* Place central faces */
    table = malloc(sizeof(int) * 8);
    for( i = 0; i < 8; ++i )
    {
        table[i] = BT_BTp + i * iedges;
    }
    faces[nfaces - (has_culet ? 2 : 1)] = make_face_n2(table, 8, 0, 0);
    if( has_culet )
    {
        int* culet = malloc(sizeof(int) * 8);
        for( i = 0; i < 8; ++i )
        {
            culet[i] = PB_C + (7-i) * iedges;
        }
        faces[nfaces - 1] = make_face_n2(culet, 8, 0, 0);
    }

    for( i = 0; i < ifaces; ++i )
    {
        free(lf[i].edges);
    }
    free(lf);

    retval = malloc( sizeof(convexpolyhedron_t) );
    *retval = convexpolyhedron_init(points, npoints, edges, nedges, faces, nfaces);

    if( !has_error )
    {
        return retval;
    }
    else
    {
        dscan_convexpolyhedron_destroy( retval );
        return 0;
    }
}

int dscan_roundcut_rotations()
{
    return 2;
}

matrix3_t dscan_roundcut_rotation
(
 int index
)
{
    if( 0 == index )
    {
        return make_yrotation(index);
    }
    else
    {
        matrix3_t r;
        double cosq = sqrt( 0.5 + sqrt(2.0) / 4.0 );
        double sinq = sqrt( 0.5 - sqrt(2.0) / 4.0 );
        assert( 1 == index );
        r.d[1] = r.d[3] = r.d[5] = r.d[7] = 0.0F;
        r.d[4] = 1.0F;
        r.d[0] = r.d[8] = cosq;
        r.d[2] = sinq;
        r.d[6] = -sinq;
        return r;
    }
}

convexpolygon_t* dscan_roundcut_crosssection
(
 convexpolyhedron_t* poly,
 int index
)
{
    bool has_culet = ((poly->mPointsCnt % 8) == 0);
    int ipoints = poly->mPointsCnt / 8;
    if( 0 == index )
    {
        if( has_culet )
        {
            vector2_t* points = malloc( sizeof(vector2_t) * 8 );
            vector3_t rpavilion_mid = mid3( poly->mPoints[PB],
                                            poly->mPoints[PB + ipoints * 7] );
            vector3_t lpavilion_mid = mid3( poly->mPoints[PB + ipoints * 4],
                                            poly->mPoints[PB + ipoints * 3] );
            copy_down( points + 0, rpavilion_mid );
            copy_down( points + 1, poly->mPoints[PT] );
            copy_down( points + 2, poly->mPoints[BB] );
            copy_down( points + 3, poly->mPoints[BT] );
            copy_down( points + 4, poly->mPoints[BT + ipoints * 4] );
            copy_down( points + 5, poly->mPoints[BB + ipoints * 4] );
            copy_down( points + 6, poly->mPoints[PT + ipoints * 4] );
            copy_down( points + 7, lpavilion_mid );
            return dscan_convexpolygon_create( points, 8 );
        }
        else
        {
            vector2_t* points = malloc( sizeof(vector2_t) * 7 );
            copy_down( points + 0, poly->mPoints[8 * ipoints] );
            copy_down( points + 1, poly->mPoints[PT] );
            copy_down( points + 2, poly->mPoints[BB] );
            copy_down( points + 3, poly->mPoints[BT] );
            copy_down( points + 4, poly->mPoints[BT + ipoints * 4] );
            copy_down( points + 5, poly->mPoints[BB + ipoints * 4] );
            copy_down( points + 6, poly->mPoints[PT + ipoints * 4] );
            return dscan_convexpolygon_create( points, 7 );
        }
    }
    else
    {
        vector2_t* points = malloc( sizeof(vector2_t) * (11+(has_culet?1:0)) );
        vector3_t rstar_mid = mid3( poly->mPoints[BT], 
                                    poly->mPoints[BT + ipoints] );
        vector3_t lstar_mid = mid3( poly->mPoints[BT + ipoints * 4],
                                    poly->mPoints[BT + ipoints * 5] );
        assert( 1 == index );
        copy_down( points + 0, poly->mPoints[has_culet? PB : 8*ipoints] );
        copy_down( points + 1, poly->mPoints[PR] );
        copy_down( points + 2, poly->mPoints[L0R] );
        copy_down( points + 3, poly->mPoints[U0R] );
        copy_down( points + 4, poly->mPoints[BR] );
        copy_down( points + 5, rstar_mid );
        copy_down( points + 6, lstar_mid );
        copy_down( points + 7, poly->mPoints[BR + ipoints * 4] );
        copy_down( points + 8, poly->mPoints[U0R + ipoints * 4] );
        copy_down( points + 9, poly->mPoints[L0R + ipoints * 4] );
        copy_down( points + 10, poly->mPoints[PR + ipoints * 4] );
        if( has_culet )
        {
            copy_down( points + 11, poly->mPoints[PB + ipoints * 4] );
        }
        return dscan_convexpolygon_create( points, 11 + (has_culet?1:0) );
    }
}
