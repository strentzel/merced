class Vector2
{
    public float x;
    public float y;

    public Vector2()
    {
        x = 0;
        y = 0;
    }

    public Vector2( Vector2 that )
    {
        x = that.x;
        y = that.y;
    }

    public Vector2( float xi, float yi )
    {
        x = xi;
        y = yi;
    }

    public float Dot( Vector2 that )
    {
        return x * that.x + y * that.y;
    }

    public float Cross( Vector2 that )
    {
        return x * that.y - y * that.x;
    }

    public float Norm()
    {
        return (float) Math.sqrt( x * x + y * y );
    }

    public float Norm2()
    {
        return x * x + y * y;
    }

    public String toString()
    {
        return "<" + x + "," + y + ">";
    }

    public float ComputeIntersectScale( Vector2 rp, Vector2 v )
    {
        //rp + a*v = b*this

        return (rp.x*v.y - v.x*rp.y) / (x*v.y - v.x*y); 
    } 

    static public void Copy( Vector2 retval, Vector2 v )
    {
        retval.x = v.x;
        retval.y = v.y;
    }

    static public void Normalize( Vector2 v )
    {
        float scale = 1.0F / (float)Math.sqrt( v.x * v.x + v.y * v.y );
        v.x *= scale;
        v.y *= scale;
    }

    static public void Scale( Vector2 result, Vector2 v, float s )
    {
        result.x = v.x * s;
        result.y = v.y * s;
    }

    static public void Add( Vector2 result, Vector2 s, Vector2 t )
    {
        result.x = s.x + t.x;
        result.y = s.y + t.y;
    }

    static public void Minus( Vector2 result, Vector2 s, Vector2 t )
    {
        result.x = s.x - t.x;
        result.y = s.y - t.y;
    }

    static public void NormFromRel( Vector2 result, Vector2 f, Vector2 r )
    {
        float fnn = f.Norm2();
        float dot = f.Dot( r );
        result.x = r.x - dot*f.x/fnn;
        result.y = r.y - dot*f.y/fnn;
        if( result.Norm2() < 1e-8F )
        {
            result.x = f.y;
            result.y = f.x;
        }
        Normalize( result );
    }
}