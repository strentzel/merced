package com.example.dscan;

public class ConvexPolyhedron
{
    ConvexPolyhedron( long impl )
    {
        mImpl = impl;
    }

    protected void finalize()
    {
        Native.ConvexPolyhedronDestroy( mImpl );
    }

    public void Transform( Matrix3 M )
    {
        Native.ConvexPolyhedronTransform( mImpl, 
                                          M.d00, M.d01, M.d02,
                                          M.d10, M.d11, M.d12,
                                          M.d20, M.d21, M.d22 );
    }

    public float ComputeVolume()
    {
        return Native.ConvexPolyhedronVolume( mImpl );
    }

    public Edge[] GetVisibleEdges()
    {
        float[] points = Native.ConvexPolyhedronVisibleEdges( mImpl );
        Edge[] retval = new Edge[points.length/4];
        for( int i = 0; i < points.length; i += 4 ){
            retval[i/4] = new Edge( new Vector2( points[i], points[i+1] ),
                                    new Vector2( points[i+2], points[i+3] ) );
        }
        return retval;
    }

    long mImpl;
}
