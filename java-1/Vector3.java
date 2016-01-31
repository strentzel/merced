class Vector3
{
    public float x;
    public float y;
    public float z;

    public Vector3()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    public Vector3( Vector3 that )
    {
        x = that.x;
        y = that.y;
        z = that.z;
    }

    public Vector3( float xi, float yi, float zi )
    {
        x = xi;
        y = yi;
        z = zi;
    }

    public float Dot( Vector3 that )
    {
        return x * that.x + y * that.y + z * that.z;
    }

    public float Cross( Vector3 that )
    {
        return y * that.z - z * that.y +
            z * that.x - x * that.z +
            x * that.y - y * that.x;
    }

    public float Norm()
    {
        return (float) Math.sqrt( x * x + y * y + z * z );
    }

    public float Norm2()
    {
        return x * x + y * y + z * z;
    }

    public String toString()
    {
        return "<" + x + "," + y + "," + z + ">";
    }

    static public void Copy( Vector3 retval, Vector3 v )
    {
        retval.x = v.x;
        retval.y = v.y;
        retval.z = v.z;
    }

    static public void Normalize( Vector3 v )
    {
        float scale = 1.0F / (float)Math.sqrt( v.x*v.x + v.y*v.y + v.z*v.z);
        v.x *= scale;
        v.y *= scale;
        v.z *= scale;
    }

    static public void Cross( Vector3 result, Vector3 s, Vector3 t )
    {
        result.x = s.y * t.z - s.z * t.y;
        result.y = s.z * t.x - s.x * t.z;
        result.z = s.x * t.y - s.y * t.x;
    }

    static public void Scale( Vector3 result, Vector3 v, float s )
    {
        result.x = v.x * s;
        result.y = v.y * s;
        result.z = v.z * s;
    }

    static public void Add( Vector3 result, Vector3 s, Vector3 t )
    {
        result.x = s.x + t.x;
        result.y = s.y + t.y;
        result.z = s.z + t.z;
    }

    static public void Minus( Vector3 result, Vector3 s, Vector3 t )
    {
        result.x = s.x - t.x;
        result.y = s.y - t.y;
        result.z = s.z - t.z;
    }

    static public void NormFromRel( Vector3 result, Vector3 f, Vector3 r )
    {
        float fnn = f.Norm2();
        float dot = f.Dot( r );
        result.x = r.x - dot*f.x/fnn;
        result.y = r.y - dot*f.y/fnn;
        result.z = r.z - dot*f.z/fnn;
        if( result.Norm2() < 1e-8F )
        {
            result.x = f.y;
            result.y = f.x;
            result.z = f.z;
        }
        Normalize( result );
    }

    static public void IntersectPlane( Vector3 result, Vector3 planeNorm,
                                       Vector3 rpoint, Vector3 rvec )
    {
        float k = planeNorm.Dot( rpoint );
        float y = planeNorm.Dot( rvec );
        Vector3.Scale( result, rvec, -k / y );
        Vector3.Add( result, result, rpoint );
    }

    static public void Normal( Vector3 result, 
                               Vector3 a, Vector3 b, Vector3 c )
    {
        float x0 = b.x - a.x;
        float y0 = b.y - a.y;
        float z0 = b.z - a.z;
        float x1 = c.x - b.x;
        float y1 = c.y - b.y;
        float z1 = c.z - b.z;
        result.x = y0*z1 - y1*z0;
        result.y = z0*x1 - z1*x0;
        result.z = x0*y1 - x1*y0;
        Normalize( result );
    }
}
