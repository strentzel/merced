package com.example.dscan;

public class Observer
{
    public Observer( String name )
    {
        mImpl = Native.ObserverCreate( name );
    }

    protected void finalize()
    {
        Native.ObserverDestroy( mImpl );
    }

    long mImpl;
}
