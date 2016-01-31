package com.example.dscan;

public class Material
{
    public Material( String name )
    {
        mImpl = Native.MaterialCreate( name );
    }

    protected void finalize()
    {
        Native.MaterialDestroy( mImpl );
    }

    long mImpl;
}

