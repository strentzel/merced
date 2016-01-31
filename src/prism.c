#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dscan.h"
#include "dscanimp.h"

static int MINWAVELENGTH = 380;
static int MAXWAVELENGTH = 780;
static int INCREMENT = 5;
static float INCREMENTINV = 0.2F; /*1.0F / (float) INCREMENT*/
static float XYZSCALE = 4.67912874e-2F; /*1.0F / 21.3715F*/
static float OBSERVER[] =
{
    /*380*/ 0.001368,0.000039,0.00645,
    /*385*/ 0.002236,0.000064,0.01055,
    /*390*/ 0.004243,0.00012,0.02005,
    /*395*/ 0.00765,0.000217,0.03621,
    /*400*/ 0.01431,0.000396,0.06785,
    /*405*/ 0.02319,0.00064,0.1102,
    /*410*/ 0.04351,0.00121,0.2074,
    /*415*/ 0.07763,0.00218,0.3713,
    /*420*/ 0.13438,0.004,0.6456,
    /*425*/ 0.21477,0.0073,1.03905,
    /*430*/ 0.2839,0.0116,1.3856,
    /*435*/ 0.3285,0.01684,1.62296,
    /*440*/ 0.34828,0.023,1.74706,
    /*445*/ 0.34806,0.0298,1.7826,
    /*450*/ 0.3362,0.038,1.77211,
    /*455*/ 0.3187,0.048,1.7441,
    /*460*/ 0.2908,0.06,1.6692,
    /*465*/ 0.2511,0.0739,1.5281,
    /*470*/ 0.19536,0.09098,1.28764,
    /*475*/ 0.1421,0.1126,1.0419,
    /*480*/ 0.09564,0.13902,0.81295,
    /*485*/ 0.05795,0.1693,0.6162,
    /*490*/ 0.03201,0.20802,0.46518,
    /*495*/ 0.0147,0.2586,0.3533,
    /*500*/ 0.0049,0.323,0.272,
    /*505*/ 0.0024,0.4073,0.2123,
    /*510*/ 0.0093,0.503,0.1582,
    /*515*/ 0.0291,0.6082,0.1117,
    /*520*/ 0.06327,0.71,0.07825,
    /*525*/ 0.1096,0.7932,0.05725,
    /*530*/ 0.1655,0.862,0.04216,
    /*535*/ 0.22575,0.91485,0.02984,
    /*540*/ 0.2904,0.954,0.0203,
    /*545*/ 0.3597,0.9803,0.0134,
    /*550*/ 0.43345,0.99495,0.00875,
    /*555*/ 0.51205,1,0.00575,
    /*560*/ 0.5945,0.995,0.0039,
    /*565*/ 0.6784,0.9786,0.00275,
    /*570*/ 0.7621,0.952,0.0021,
    /*575*/ 0.8425,0.9154,0.0018,
    /*580*/ 0.9163,0.87,0.00165,
    /*585*/ 0.9786,0.8163,0.0014,
    /*590*/ 1.0263,0.757,0.0011,
    /*595*/ 1.0567,0.6949,0.001,
    /*600*/ 1.0622,0.631,0.0008,
    /*605*/ 1.0456,0.5668,0.0006,
    /*610*/ 1.0026,0.503,0.00034,
    /*615*/ 0.9384,0.4412,0.00024,
    /*620*/ 0.85445,0.381,0.00019,
    /*625*/ 0.7514,0.321,0.0001,
    /*630*/ 0.6424,0.265,0.00005,
    /*635*/ 0.5419,0.217,0.00003,
    /*640*/ 0.4479,0.175,0.00002,
    /*645*/ 0.3608,0.1382,0.00001,
    /*650*/ 0.2835,0.107,0,
    /*655*/ 0.2187,0.0816,0,
    /*660*/ 0.1649,0.061,0,
    /*665*/ 0.1212,0.04458,0,
    /*670*/ 0.0874,0.032,0,
    /*675*/ 0.0636,0.0232,0,
    /*680*/ 0.04677,0.017,0,
    /*685*/ 0.0329,0.01192,0,
    /*690*/ 0.0227,0.00821,0,
    /*695*/ 0.01584,0.005723,0,
    /*700*/ 0.011359,0.004102,0,
    /*705*/ 0.008111,0.002929,0,
    /*710*/ 0.00579,0.002091,0,
    /*715*/ 0.004109,0.001484,0,
    /*720*/ 0.002899,0.001047,0,
    /*725*/ 0.002049,0.00074,0,
    /*730*/ 0.00144,0.00052,0,
    /*735*/ 0.001,0.000361,0,
    /*740*/ 0.00069,0.000249,0,
    /*745*/ 0.000476,0.000172,0,
    /*750*/ 0.000332,0.00012,0,
    /*755*/ 0.000235,0.000085,0,
    /*760*/ 0.000166,0.00006,0,
    /*765*/ 0.000117,0.000042,0,
    /*770*/ 0.000083,0.00003,0,
    /*775*/ 0.000059,0.000021,0,
    /*780*/ 0.000042,0.000015,0
};

static int sign( float f )
{
    if( f > 0 ) return 1;
    else return -1;
}

static void material_init_air
(
 material_t* m
)
{
    m->n_min = 1.0F;
    m->n_max = 1.0F;
    m->theta = (float) M_PI / 2.0F;
    m->cos = 0.0F;
    m->c0 = m->c1 = m->c2 = 0.0F;
}

static void material_init_curve
(
 material_t* m,
 float a,
 float b,
 float c
)
{
    m->c0 = -b / (2.0F * a);
    m->c1 = 1.0F / a;
    m->c2 = (b / (2.0F * a)) * (b / (2.0F * a)) - c / a;
    m->n_min = -m->c2 / m->c1;
    m->n_max = (float) (a * 380 * 380 + b * 380 + c);
    m->theta = (float) asin(1.0 / m->n_min);
    m->cos = (float) cos(m->theta);
}

static void material_init_constant
(
 material_t* m,
 float n
)
{
    m->n_min = n - 0.01F;
    m->n_max = n + 0.01F;
    m->c0 = 731;
    m->c1 = 3.125E6;
    m->c2 = -m->c1 * m->n_min;
    m->theta = (float) asin(1.0 / m->n_min);
    m->cos = (float) cos(m->theta);
}

material_t* dscan_material_create
(
 const char* type,
 int len
)
{
    material_t* retval = malloc( sizeof(material_t) );
    if( 7 == len && 0 == strncmp( "diamond", type, len ) )
    {
        material_init_curve( retval, 4.98E-7F, -7.28E-4F, 2.67F );
    }
    else if( 11 == len && 0 == strncmp( "crownimpure", type, len ) )
    {
        material_init_constant( retval, sqrt(3.0F) );
    }
    else
    {
        material_init_air( retval );
    }
    return retval;
}

void dscan_material_destroy
(
 material_t* mat
)
{
    free( mat );
}

observer_t* dscan_observer_create
(
 const char* type,
 int len
)
{
    observer_t* retval = malloc( sizeof(observer_t) );
    retval->type = 0;
    return retval;
}

void dscan_observer_destroy
(
 observer_t* obs
)
{
    free( obs );
}

typedef vector3_t Xyz;

#define GETINDEX(w) ((((int)w)-MINWAVELENGTH) / INCREMENT)
/*
static int GetIndex( float wavelength )
{
    int iwave = (int) wavelength;
    assert( iwave >= MINWAVELENGTH && iwave <= MAXWAVELENGTH );
    return (iwave - MINWAVELENGTH) / INCREMENT;
}
*/

static vector3_t Integrate( observer_t* obs, float w0, float w1 )
{
    float* data = OBSERVER;
    vector3_t datum;
    float s;
    int i;

    assert( 0 == obs->type );

    datum.x = datum.y = datum.z = 0.0F;
    if( w0 < MAXWAVELENGTH && w1 > MINWAVELENGTH )
    {
        int i0, i1;

        if( w0 < MINWAVELENGTH )
        {
            w0 = MINWAVELENGTH;
        }
        if (w1 > MAXWAVELENGTH )
        {
            w1 = MAXWAVELENGTH;
        }

        i0 = GETINDEX(w0);
        i1 = GETINDEX(w1);

        if( i0 == i1 )
        {
            float* d = data + i0*3;
            s = (w1 - w0) * INCREMENTINV;
            datum.x = *d++ * s;
            datum.y = *d++ * s;
            datum.z = *d++ * s;
        }
        else
        {
            float* d = data + i0*3;
            s = ((MINWAVELENGTH + (i0 + 1 ) * INCREMENT) - w0) * INCREMENTINV;
            datum.x = *d++ * s;
            datum.y = *d++ * s;
            datum.z = *d++ * s;
            
            for( i = i0 + 1; i < i1; ++i ) 
            {
                datum.x += *d++;
                datum.y += *d++;
                datum.z += *d++;
            }
            
            s = (w1 - (MINWAVELENGTH + (i1) * INCREMENT)) * INCREMENTINV;
            datum.x += *d++ * s;
            datum.y += *d++ * s;
            datum.z += *d++ * s;
        }
    } 
    else
    {
        datum.x = 0.0F;
        datum.y = 0.0F;
        datum.z = 0.0F;
    }

    datum.x *= XYZSCALE;
    datum.y *= XYZSCALE;
    datum.z *= XYZSCALE;

    /*assert( datum.x >= 0 && datum.y >= 0 && datum.z >= 0 );*/
    /*datum.x = Bound(datum.x);
    datum.y = Bound(datum.y);
    datum.z = Bound(datum.z);*/

    return datum;
}

#define Ninety (M_PI / 2.0F)

static vector2_t NormFromRel( vector2_t f, vector2_t r )
{
    vector2_t result;
    float fnn = norm2( f );
    float dot = dot2( f, r );
    result.x = r.x - dot * f.x / fnn;
    result.y = r.y - dot * f.y / fnn;
    if( norm2( result ) < 1e-8F )
    {
        result.x = f.y;
        result.y = f.x;
    }
    return normalize2( result );
}

static int cap(float f) {
    if (f < 0)
        return 0;
    else if (f > 254.5F)
        return 255;
    else
        return (int) (f + 0.5F);
}

static Xyz toXyz( int argb )
{
    float a = ((argb >> 24) & 0xFF) / 255.0F;
    int r = ((argb >> 16) & 0xFF);
    int g = ((argb >> 8) & 0xFF);
    int b = ((argb) & 0xFF);
    Xyz retval;
    float m;
    if( 0 == r && 0 == g && 0 == b )
    {
        m = 0;
        a = 0;
    }
    else
    {
        assert( r > 0 || g > 0 || b > 0 );
        assert( r >= 0 && g >= 0 && b >= 0 );
        m = a / (0.17697 * r + 0.8124 * g + 0.01063 * b);
    }
    retval.x = (0.49F * r + 0.31 * g + 0.2 * b) * m;
    retval.y = a;
    retval.z = (0.01 * g + 0.99 * b) * m;
    return retval;
}

static int toArgb( Xyz xyz, float weight )
{
    float rf = 2.364614F * xyz.x - 0.896541F * xyz.y - 0.468073F * xyz.z;
    float gf = -0.515166F * xyz.x + 1.426408F * xyz.y + 0.088758F * xyz.z;
    float bf = 0.005204F * xyz.x - 0.014408F * xyz.y + 1.009204F * xyz.z;
    float mx = rf;
    float af;
    int a, r, g, b;
    assert( weight >= 0.0F && weight <= 1.0F );

    if (gf > mx)
        mx = gf;
    if (bf > mx)
        mx = bf;
    if( mx <= 0 ) return 0;
    af = (xyz.y*weight < 1.0F ? xyz.y*weight : 1.0F);
    a = cap(af * 255.0F);
    r = cap(rf / mx * 255.0F);
    g = cap(gf / mx * 255.0F);
    b = cap(bf / mx * 255.0F);
    return (a << 24) | (r << 16) | (g << 8) | b;
}

#define NORM2(v) sqrt(v.x*v.x+v.y*v.y)
#define DOT2(a,b) (a.x*b.x+a.y*b.y)

static int ComputeFull( prism_t* prism, vector2_t point )
{
    float first, last;
    int intvls = 10;
    float dist, weight, mult;
    Xyz xyz;
    int i;
    vector2_t tmp1, tmp2;
    bool err;
/*    vector2_t tmp1 = minus2( point, prism->mSurfacePoint );
    vector2_t tmp2 = scale2( prism->mNormB, -0.5F );
    tmp1 = add2( tmp1, tmp2 );
    first = (float) min(1.0, max(COMPUTEINTERSECTSCALE( prism->mSurfaceVec, tmp1, prism->mNegB), 0.0));*/
    tmp1.x = (point.x - prism->mSurfacePoint.x) + prism->mNormB.x * -0.5F;
    tmp1.y = (point.y - prism->mSurfacePoint.y) + prism->mNormB.y * -0.5F;
/*    num = tmp1.x * prism->mNegB.y - prism->mNegB.x * tmp1.y;
    den = prism->mSurfaceVec.x * prism->mNegB.y + prism->mNegB.x * prism->mSurfaceVec.y;
    if( den < 0 ){ den *= -1; num *= -1; }
    assert( den >= 0 );
    if( num > den ) first = 1;
    else if( num < 0 ) first = 0;
    else first = num / den;*/
    first = compute_intersect_scale( prism->mSurfaceVec, tmp1, prism->mNegB,
                                     &err );
/*    tmp1 = minus2( point, prism->mSurfacePoint );
    tmp2 = scale2( prism->mNormC, -0.5F );
    tmp1 = add2( tmp1, tmp2 );*/
    tmp1.x = (point.x - prism->mSurfacePoint.x) + prism->mNormC.x * -0.5F;
    tmp1.y = (point.y - prism->mSurfacePoint.y) + prism->mNormC.y * -0.5F;
/*    last = (float) max(0.0, min(COMPUTEINTERSECTSCALE( prism->mSurfaceVec, tmp1, prism->mNegC), 1.0));*/
    last = compute_intersect_scale( prism->mSurfaceVec, tmp1, prism->mNegC,
                                    &err );

    dist = (last - first) * NORM2( prism->mSurfaceVec );
    weight = dist / intvls;
    mult = (last - first) / intvls;
    xyz.x = 0.0F;
    xyz.y = 0.0F;
    xyz.z = 0.0F;
    /* TODO Large distances should reduce need to sub-divide */
    for( i = 0; i < intvls; ++i )
    {
        float distance, cosTheta, sinTheta0, sinTheta1;
        float m = first + (i + 0.5F) * mult;
        float x, pre, post;
        /*tmp1 = scale2( prism->mSurfaceVec, first + (i + 0.5F) * mult);
        tmp1 = add2( prism->mSurfacePoint, tmp1 );
        tmp2 = minus2( point, tmp1 ); */
        tmp2.x = point.x - (prism->mSurfacePoint.x + prism->mSurfaceVec.x*m);
        tmp2.y = point.y - (prism->mSurfacePoint.y + prism->mSurfaceVec.y*m);
        distance = NORM2( tmp2 );
        /*theta = (float) acos( DOT2( tmp2, prism->mNormE ) / distance);
        theta0 = theta - (float) asin(0.5 / distance);
        theta1 = theta + (float) asin(0.5 / distance);*/
        cosTheta = DOT2( tmp2, prism->mNormE ) / distance;
        x = 0.5 / distance;
        pre = sqrt( (1-cosTheta*cosTheta)*(1 - x*x) );
        post = cosTheta * x;
        sinTheta0 = pre - post;
        sinTheta1 = pre + post;

        if((sinTheta0 < prism->mSinMaxTheta || sinTheta1 > prism->mSinMinTheta) && DOT2( tmp2, prism->mSurfaceVec) > 0)
        {
            float n0 = (float) max(prism->material_incident->n_min, (sinTheta0 / prism->mSinThetaD));
            float n1 = (float) min(prism->material_incident->n_max, (sinTheta1 / prism->mSinThetaD));

            float c0 = prism->material_incident->c0;
            float c1 = prism->material_incident->c1;
            float c2 = prism->material_incident->c2;

            float w0 = (float) (c0 - sqrt(max(c1 * n0 + c2, 0)));
            float w1 = (float) (c0 - sqrt(max(c1 * n1 + c2, 0)));

            float n = (n0 + n1) / 2;
            float nx = (prism->material_transmittent->n_min +
                        prism->material_transmittent->n_max) / 2;
            float A = 2 * n * nx * prism->mCosThetaD * cosTheta;
            float ds = n * prism->mCosThetaD + nx * cosTheta;
            float dp = n * cosTheta + nx * prism->mCosThetaD;
            float T = A / (ds*ds) + A / (dp*dp);
            float maxTd = 2 * n / nx + 2 * nx / n;
            float Td;
            float adj_weight;
            Xyz tmp;

            if( T * prism->mCosThetaD > maxTd * cosTheta )
            {
                Td = maxTd;
            }
            else
            {
                Td = T * prism->mCosThetaD / cosTheta;
            }
            adj_weight = prism->light_weight * weight * Td;

            /* Check on bounding in Integrate.  Why negative? */
            tmp = Integrate( prism->observer, w1, w0);
            xyz.x += tmp.x * adj_weight;
            xyz.y += tmp.y * adj_weight;
            xyz.z += tmp.z * adj_weight;
        }
    }

    return toArgb(xyz, 1.0F);
}

static int GetRegion( prism_t* prism, vector2_t point )
{
    /*vector2_t tmp1 = minus2( point, prism->mSurfacePoint );*/
    vector2_t tmp1;
    tmp1.x = point.x - prism->mSurfacePoint.x;
    tmp1.y = point.y - prism->mSurfacePoint.y;
    if( DOT2( tmp1, prism->mNormE ) > 1E-4F )
    {
        /*vector2_t tmp2 = minus2( tmp1, prism->mSurfaceVec );*/
        vector2_t tmp2;
        tmp2.x = tmp1.x - prism->mSurfaceVec.x;
        tmp2.y = tmp1.y - prism->mSurfaceVec.y;

        if( DOT2( tmp1, prism->mNormB ) >= 1 && DOT2( tmp2, prism->mNormC ) >= 1 )
        {
            assert( DOT2( tmp1, prism->mNormA ) >= -1 &&
                    DOT2( tmp2, prism->mNormD ) >= -1 );
            return 3;
        }
        else if( DOT2( tmp1, prism->mNormA ) < -1 || DOT2( tmp2, prism->mNormD) < -1)
        {
            return 2;
        }
        else
        {
            return 4;
        }
    }
    else
    {
        return 1;
    }
}

static int Compute( prism_t* prism, vector2_t point, int region) {
    switch (region) {
    case 1:
    case 2:
        return 0x00000000;
    case 3:
        if( prism->mHasCoreSample )
        {
            return prism->mCoreRgb;
        }
        else
        {
            return ComputeFull( prism, point );
        }
    case 4:
        return ComputeFull( prism, point );
    }
    return 0x00000000;
}

static int GetPixel( prism_t* prism, int x, int y )
{
    vector2_t point;
    point.x = x;
    point.y = y;
    return Compute( prism, point, GetRegion(prism, point));
}

static convexpolygon_t IntersectPrism( prism_t* prism, convexpolygon_t* bounds )
{
    vector2_t tmp1, tmp2, tmp3;
    vector2_t* points;
    float radius, length;
    convexpolygon_t poly, retval;
    int cnt;

    if( prism->mSinMinTheta >= 1.0F )
    {
        assert( prism->mSinMaxTheta >= 1.0F );
        return convexpolygon_init( 0, 0 );
    }

    tmp1 = scale2( prism->mSurfaceVec, 0.5F );
    tmp1 = add2( prism->mSurfacePoint, tmp1 );

    tmp2 = bounds->center;
    tmp2 = minus2( tmp2, tmp1);
    radius = sqrt( bounds->radius2 ) + norm2( tmp2 );

    tmp1 = scale2( prism->mSurfaceVec, 0.5F );
    tmp1 = minus2( tmp1, prism->mMinDir);

    tmp2 = add2( prism->mSurfaceVec, prism->mMaxDir );
    tmp2 = normalize2( minus2( tmp2, prism->mMinDir ) );

    tmp2 = scale2( tmp2, dot2( tmp1, tmp2) );
    tmp1 = minus2( tmp1, tmp2 );

    length = radius / norm2( tmp1 );
    length = (float) min(length, 1000);

    tmp1 = add2( prism->mSurfacePoint, prism->mSurfaceVec );
    tmp2 = scale2( prism->mMaxDir, length );
    tmp2 = add2( tmp1, tmp2 );
    tmp3 = scale2( prism->mMinDir, length );
    tmp3 = add2( prism->mSurfacePoint, tmp3 );

    if( prism->mSinMaxTheta <= 1.0F )
    {
        points = malloc( sizeof(vector2_t) * 4 );
        cnt = 4;
        if( cross2( prism->mSurfaceVec, prism->mNormE ) > 0 )
        {
            points[0] = prism->mSurfacePoint;
            points[1] = tmp1;
            points[2] = tmp2;
            points[3] = tmp3;
        }
        else
        {
            points[3] = prism->mSurfacePoint;
            points[2] = tmp1;
            points[1] = tmp2;
            points[0] = tmp3;
        }
    }
    else
    {
        points = malloc( sizeof(vector2_t) * 3 );
        cnt = 3;
        if( cross2( prism->mSurfaceVec, prism->mNormE ) > 0 )
        {
            points[0] = prism->mSurfacePoint;
            points[1] = tmp2;
            points[2] = tmp3;
        }
        else
        {
            points[2] = prism->mSurfacePoint;
            points[1] = tmp2;
            points[0] = tmp3;
        }
    }

    poly = convexpolygon_init( points, cnt );
    retval = convexpolygon_intersect( &poly, bounds );
    convexpolygon_destroy( poly );
    return retval;
}

static void addToPixel( int* buf, int argb )
{
/*
    int curr = *buf;
    int adst = ((curr >> 24) & 0xFF);
    int asrc = ((argb >> 24) & 0xFF);
    int ap =  adst + asrc;
    int rp, gp, bp;
    if( ap <= 0 ) return;
    rp = (((curr >> 16) & 0xFF)*adst + ((argb >> 16) & 0xFF)*asrc)/ap;
    gp = (((curr >> 8) & 0xFF)*adst + ((argb >> 8) & 0xFF)*asrc)/ap;
    bp = (((curr) & 0xFF)*adst + ((argb) & 0xFF)*asrc)/ap;
    if( ap > 255 ) ap = 255;
    if( rp > 255 ) rp = 255;
    if( gp > 255 ) gp = 255;
    if( bp > 255 ) bp = 255;
    *buf = (ap << 24) | (rp << 16) | (gp << 8) | bp;
    */
    Xyz dst = toXyz( *buf );
    Xyz src = toXyz( argb );
    *buf = toArgb( add3( dst, src ), 1.0F );
}

static convexpolygon_t ComputeBounds( prism_t* prism, int w, int h )
{
    vector2_t* points;
    convexpolygon_t screen;

    points = malloc( sizeof(vector2_t) * 4 );
    points[0].x = 0.0;
    points[0].y = 0.0;
    points[1].x = w;
    points[1].y = 0.0;
    points[2].x = w;
    points[2].y = h;
    points[3].x = 0.0;
    points[3].y = h;
    screen = convexpolygon_init( points, 4 );
    return IntersectPrism( prism, &screen );
}

prism_t* dscan_prism_create( vector2_t surface1_raw, vector2_t surface2_raw, vector2_t lightDir_raw, material_t* material_incident, material_t* material_transmittent, observer_t* observer, matrix3_t display_transform, float weight )
{
    prism_t* retval = malloc( sizeof(prism_t) );
    vector2_t tmp1;
    vector2_t tmp2;
    vector2_t surfaceDir, normal;
    float coreSampleLen;
    bool extreme;
    vector2_t surface1, surface2, lightDir;
    bool err;
    surface1 = multiply3( display_transform, surface1_raw, 1 );
    surface2 = multiply3( display_transform, surface2_raw, 1 );
    lightDir = multiply3( display_transform, lightDir_raw, 0 );

    retval->material_incident = material_incident;
    retval->material_transmittent = material_transmittent;
    retval->observer = observer;
    retval->light_weight = weight;

    tmp1 = minus2( surface2, surface1 );
    if( dot2( lightDir, tmp1 ) >= 0.0F )
    {
        retval->mSurfacePoint = surface1;
        tmp1 = minus2( surface2, surface1 );
        retval->mSurfaceVec = tmp1;
    }
    else
    {
        retval->mSurfacePoint = surface2;
        tmp1 = minus2( surface1, surface2 );
        retval->mSurfaceVec = tmp1;
    }

    surfaceDir = normalize2( retval->mSurfaceVec );
    retval->mThetaD = Ninety - (float) acos( dot2( surfaceDir, lightDir ) / norm2( lightDir ));
    normal.x =  surfaceDir.y;
    normal.y = -surfaceDir.x;
    normalize2( normal );
    if( dot2( lightDir, normal ) < 0 )
    {
        normal = scale2( normal, -1 );
    }

    retval->mNormE = NormFromRel( surfaceDir, lightDir );
    retval->mSinMinTheta = material_incident->n_min * sin(retval->mThetaD);

    extreme = 1;
    if( retval->mSinMinTheta < 1.0F )
    {
        tmp1 = scale2( surfaceDir, retval->mSinMinTheta );
        tmp2 = scale2( normal, sqrt( 1.0 - retval->mSinMinTheta*retval->mSinMinTheta) );
        retval->mMinDir = normalize2( add2( tmp1, tmp2 ) );

        tmp1 = scale2( surfaceDir, -1 );
        retval->mNormA = NormFromRel( retval->mMinDir, surfaceDir );
        retval->mNormC = NormFromRel( retval->mMinDir, tmp1 );

        /* make sure that errors do not cause minDir to hit the incident side*/
        if( sign( cross2( normal, surfaceDir ) ) ==
            sign( cross2( retval->mMinDir, surfaceDir ) ) )
        {
            extreme = 0;
        }
    }
    if( extreme )
    {
        retval->mSinMinTheta = 1.0F;
        retval->mMinDir = surfaceDir;
        tmp1 = scale2( lightDir, -1 );
        retval->mNormA = NormFromRel( retval->mMinDir, tmp1 );
        retval->mNormC = NormFromRel( retval->mMinDir, lightDir );
    }
    retval->mNegC = scale2( retval->mMinDir, -1 / norm2( retval->mMinDir) );

    retval->mSinMaxTheta = material_incident->n_max * sin(retval->mThetaD);
    extreme = 1;
    if( retval->mSinMaxTheta < 1.0F )
    {
        tmp1 = scale2( surfaceDir, retval->mSinMaxTheta );
        tmp2 = scale2( normal, sqrt( 1.0F - retval->mSinMaxTheta * retval->mSinMaxTheta) );
        retval->mMaxDir = normalize2( add2( tmp1, tmp2 ) );

        tmp1 = scale2( surfaceDir, -1 );
        retval->mNormB = NormFromRel( retval->mMaxDir, surfaceDir );
        retval->mNormD = NormFromRel( retval->mMaxDir, tmp1 );

        /* make sure that errors do not cause minDir to hit the incident side*/
        if( sign( cross2( normal, surfaceDir ) ) ==
            sign( cross2( retval->mMaxDir, surfaceDir ) ) )
        {
            extreme = 0;
        }
    }
    if( extreme )
    {
        retval->mSinMaxTheta = 1.0F;
        retval->mMaxDir = surfaceDir;

        tmp1 = scale2( lightDir, -1 );
        retval->mNormB = NormFromRel( retval->mMaxDir, tmp1 );
        retval->mNormD = NormFromRel( retval->mMaxDir, lightDir );
    }
    retval->mNegB = scale2( retval->mMaxDir, -1 / norm2( retval->mMaxDir ) );

    retval->mSinThetaD = sin( retval->mThetaD );
    retval->mCosThetaD = cos( retval->mThetaD );

    coreSampleLen = compute_intersect_scale( retval->mMaxDir,
                                             retval->mSurfaceVec, 
                                             scale2(retval->mMinDir,50.0F),
                                             &err ) * 25.0F;
/*    region3Tip = COMPUTEINTERSECTSCALE( retval->mMaxDir, retval->mSurfaceVec, retval->mMinDir );
    if (fabs(region3Tip) < 50) {
        coreSampleLen = (float) fabs(region3Tip) / 2.0F;
    } else {
        coreSampleLen = 50;
        }*/

    tmp1 = scale2( retval->mSurfaceVec, 0.5F );
    tmp1 = add2( tmp1, retval->mSurfacePoint );
    tmp2 = add2( retval->mMaxDir, retval->mMinDir );
    tmp2 = scale2( tmp2, coreSampleLen );
    tmp1 = add2( tmp1, tmp2);
    if( 3 == GetRegion(retval, tmp1) )
    {
        retval->mHasCoreSample = true;
        retval->mCoreRgb = ComputeFull( retval, tmp1 );
    }
    else
    {
        retval->mHasCoreSample = false;
        retval->mCoreRgb = 0x00000000;
    }

    return retval;
}

void dscan_prism_destroy( prism_t* prism )
{
    free( prism );
}

void dscan_prism_rasterize( prism_t* prism,
                            int* bitmap, int offset, int stride, int w, int h )
{
    int lowestLeftIndex = 0;
    int highestLeftIndex = 0;
    int lowestRightIndex = 0;
    int highestRightIndex = 0;
    int i, y;
    int yStart, yEnd;
    int left, nextLeft, right, nextRight;
    int nextMinX, nextMaxX;
    convexpolygon_t poly = ComputeBounds( prism, w, h );

    if( poly.mPointsCnt <= 0 )
    {
        convexpolygon_destroy( poly );
        return;
    }

    for( i = 0; i < poly.mPointsCnt; ++i )
    {
        float x = poly.mPoints[i].x;
        float y = poly.mPoints[i].y;
        float llx = poly.mPoints[lowestLeftIndex].x;
        float lly = poly.mPoints[lowestLeftIndex].y;
        float hlx = poly.mPoints[highestLeftIndex].x;
        float hly = poly.mPoints[highestLeftIndex].y;
        float lrx = poly.mPoints[lowestRightIndex].x;
        float lry = poly.mPoints[lowestRightIndex].y;
        float hrx = poly.mPoints[highestRightIndex].x;
        float hry = poly.mPoints[highestRightIndex].y;
        if( y < lly || (y == lly && x < llx) ) lowestLeftIndex = i;
        if( y > hly || (y == hly && x < hlx) ) highestLeftIndex = i;
        if( y < lry || (y == lry && x > lrx) ) lowestRightIndex = i;
        if( y > hry || (y == hry && x > hrx) ) highestRightIndex = i;
    }

    yStart = (int) poly.mPoints[lowestLeftIndex].y;
    yEnd = (int) poly.mPoints[highestRightIndex].y;

    left = lowestLeftIndex;
    nextLeft = (left + poly.mPointsCnt - 1) % poly.mPointsCnt;
    right = lowestRightIndex;
    nextRight = (right + 1) % poly.mPointsCnt;


    nextMinX = poly.mPoints[left].x;
    nextMaxX = poly.mPoints[right].x;

    if( yEnd >= h ) yEnd = h-1;

    for( y = yStart; y <= yEnd; ++y )
    {
        int currTop = y + 1;
        int minX = nextMinX;
        int maxX = nextMaxX;

        while( poly.mPoints[nextLeft].y < currTop &&
               nextLeft != highestLeftIndex)
        {
            int nextX = (int) poly.mPoints[nextLeft].x;
            if( nextX < minX ) minX = nextX;
            left = nextLeft;
            nextLeft = (nextLeft + poly.mPointsCnt - 1) % poly.mPointsCnt;
        }

        if( y != yEnd )
        {
            float y0 = currTop;
            float sx = poly.mPoints[left].x;
            float sy = poly.mPoints[left].y;
            float tx = poly.mPoints[nextLeft].x;
            float ty = poly.mPoints[nextLeft].y;
            /* Note: computation should be valid even for tiny ty-sy since 
             * that implies tiny y0-sy since y0-sy < ty-sy */
            assert( fabs(y0-sy) <= fabs(ty-sy) );
            nextMinX = ((y0 - sy) / (ty - sy)) * (tx - sx) + sx;
        }
                                                                              
        if( nextMinX < minX ) minX = nextMinX;

        while( poly.mPoints[nextRight].y < currTop && 
               nextRight != highestRightIndex)
        {
            int nextX = (int) poly.mPoints[nextRight].x;
            if( nextX > maxX ) maxX = nextX;
            right = nextRight;
            nextRight = (nextRight + 1) % poly.mPointsCnt;
        }

        if( y != yEnd )
        {
            float y0 = currTop;
            float sx = poly.mPoints[right].x;
            float sy = poly.mPoints[right].y;
            float tx = poly.mPoints[nextRight].x;
            float ty = poly.mPoints[nextRight].y;
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
                int pixel;
                int* buf;
                assert( x >= 0 && x < w && y >= 0 && y < h );
                pixel = GetPixel( prism, x+0.5F, y+0.5F );
                buf = bitmap + offset + y * stride + x;
                addToPixel( buf, pixel );
            }
        }
    }

    convexpolygon_destroy( poly );
}
