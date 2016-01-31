#include <assert.h>
#include "dscan.h"
#include "dscanimp.h"

vector3_t normalize3( vector3_t v )
{
    float scale = 1.0F / sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
    vector3_t retval;
    retval.x = v.x * scale;
    retval.y = v.y * scale;
    retval.z = v.z * scale;
    return retval;
}

vector3_t add3( vector3_t a, vector3_t b )
{
    vector3_t retval;
    retval.x = a.x + b.x;
    retval.y = a.y + b.y;
    retval.z = a.z + b.z;
    return retval;
}

vector3_t minus3( vector3_t a, vector3_t b )
{
    vector3_t retval;
    retval.x = a.x - b.x;
    retval.y = a.y - b.y;
    retval.z = a.z - b.z;
    return retval;
}

vector3_t scale3( vector3_t v, float s )
{
    vector3_t retval;
    retval.x = v.x * s;
    retval.y = v.y * s;
    retval.z = v.z * s;
    return retval;
}

vector3_t cross3(vector3_t a, vector3_t b)
{
    vector3_t retval;
    retval.x = a.y * b.z - a.z * b.y;
    retval.y = a.z * b.x - a.x * b.z;
    retval.z = a.x * b.y - a.y * b.x;
    return retval;
}

vector3_t intersect_plane( vector3_t norm, vector3_t rpoint, vector3_t rvec )
{
    float k = dot3( norm, rpoint );
    float y = dot3( norm, rvec );
    vector3_t result = scale3( rvec, -k / y );
    result = add3( result, rpoint );
    return result;
}

float dot3( vector3_t a, vector3_t b )
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
} 

float norm3( vector3_t v ) {
    return sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
}

vector2_t multiply3( matrix3_t m, vector2_t v, bool z )
{
    vector2_t retval;
    retval.x = m.d[0] * v.x + m.d[1] * v.y;
    retval.y = m.d[3] * v.x + m.d[4] * v.y;
    if( z )
    {
        retval.x += m.d[2];
        retval.y += m.d[5];
    }
    return retval;
}

vector2_t normalize2( vector2_t v )
{
    float scale = 1.0F / sqrt( v.x*v.x + v.y*v.y );
    vector2_t retval;
    retval.x = v.x * scale;
    retval.y = v.y * scale;
    return retval;
}

vector2_t add2( vector2_t a, vector2_t b )
{
    vector2_t retval;
    retval.x = a.x + b.x;
    retval.y = a.y + b.y;
    return retval;
}

vector2_t minus2( vector2_t a, vector2_t b )
{
    vector2_t retval;
    retval.x = a.x - b.x;
    retval.y = a.y - b.y;
    return retval;
}

vector2_t scale2( vector2_t v, float s )
{
    vector2_t retval;
    retval.x = v.x * s;
    retval.y = v.y * s;
    return retval;
}

float cross2( vector2_t a, vector2_t b )
{
    return a.x*b.y - a.y*b.x;
}

float dot2( vector2_t a, vector2_t b )
{
    return a.x*b.x + a.y*b.y;
}

float norm2( vector2_t a )
{
    return sqrt( a.x*a.x + a.y*a.y );
}

float compute_intersect_scale( vector2_t s, vector2_t rp, vector2_t v,
                               bool* err ) {
    /* rp + a*v = b*s, where b is the return value.*/
    float num = rp.x * v.y - v.x * rp.y;
    float den = s.x * v.y - v.x * s.y;

    if( den > 0 && num >= den ) return 1.0F;
    else if( den < 0 && num <= den ) return 1.0F;
    else if( den > 0 && num < 0 ) return 0.0F;
    else if( den < 0 && num > 0 ) return 0.0F;
    else if( den == 0 ){ *err = 1; return 1.0F; }
    else
    {
        float retval = num / den;
        assert( retval >= 0 && retval <= 1 );
        return retval;
    }
}
