package com.example.dscan;

public class ConvexPolygon
{
    ConvexPolygon( long impl )
    {
        mImpl = impl;
    }

    protected void finalize()
    {
        Native.ConvexPolygonDestroy( mImpl );
    }

    long mImpl;
}
