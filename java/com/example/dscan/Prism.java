package com.example.dscan;

public class Prism
{
    Prism( Material external, Material internal, Observer observer,
           Vector2 start, Vector2 end, Vector2 incident, float weight )
    {
        mStart = start;
        mEnd = end;
        mIncident = incident;
        mWeight = weight;
        mExternal = external;
        mInternal = internal;
        mObserver = observer;
    }

    void Rasterize( Matrix3 displayTransform, float weight, int[] bitmap, int offset, int stride, int w, int h )
    {
        long p = Native.PrismCreate( mStart.x, mStart.y, mEnd.x, mEnd.y,
                                     mIncident.x, mIncident.y, mInternal.mImpl,
                                     mExternal.mImpl, mObserver.mImpl, 
                                     displayTransform.d00,
                                     displayTransform.d01,
                                     displayTransform.d02,
                                     displayTransform.d10,
                                     displayTransform.d11,
                                     displayTransform.d12,
                                     mWeight * weight );
        Native.PrismRasterize( p, bitmap, offset, stride, w, h );
        Native.PrismDestroy( p );
    }

    private Vector2 mStart;
    private Vector2 mEnd;
    private Vector2 mIncident;
    private float mWeight;
    private Material mExternal;
    private Material mInternal;
    private Observer mObserver;
}
