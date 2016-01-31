class Triangle
{
    Vector2 a;
    Vector2 b;
    Vector2 c;

    public Triangle( Vector2 a0, Vector2 b0, Vector2 c0 )
    {
        if( a0.x < b0.x )
        {
            if( b0.x < c0.x )
            {
                Init( a0, b0, c0 );
            }
            else if( a0.x < c0.x )
            {
                Init( a0, c0, b0 );
            }
            else if( c0.x < a0.x )
            {
                Init( c0, a0, b0 );
            }
            else if( a0.y < c0.y )
            {
                Init( a0, c0, b0 );
            }
            else
            {
                Init( c0, a0, b0 );
            }
        }
        else if( b0.x < a0.x )
        {
            if( a0.x < c0.x )
            {
                Init( b0, a0, c0 );
            }
            else if( b0.x < c0.x )
            {
                Init( b0, c0, a0 );
            }
            else if( c0.x < b0.x )
            {
                Init( c0, b0, a0 );
            }
            else if( b0.y < c0.y )
            {
                Init( b0, c0, a0 );
            }
            else
            {
                Init( c0, b0, a0 );
            }
        }
        else 
        {
            if( a0.x < c0.x )
            {
                if( a0.y < b0.y )
                {
                    Init( a0, b0, c0 );
                }
                else
                {
                    Init( b0, a0, c0 );
                }
            }
            else
            {
                Init( c0, a0, b0 );
            }
        }
    }

    private void Init( Vector2 a0, Vector2 b0, Vector2 c0 )
    {
        a = new Vector2( a0 );
        b = new Vector2( b0 );
        c = new Vector2( c0 );
    }

    public void Rasterize( RawBitmap bitmap, int rgb )
    {
        float y = (b.x-a.x)*(c.y-a.y)/(c.x-a.x) + a.y;
        Integer reflect = new Integer( (int)Math.round( c.x ) );
        float r = (float)reflect.intValue();
        Vector2 a2 = new Vector2( r + (r-c.x), c.y );
        if( b.y < y )
        {
            Vector2 b2 = new Vector2( r + (r-b.x), b.y ); 

            assert b.x >= a.x;
            if( b.x - a.x > 1.0E-3F )
            {
                Rasterize( bitmap, a, b, y, rgb, null );
            }
            if( b2.x - a2.x > 1.0E-3F )
            {
                Rasterize( bitmap, a2, b2, y, rgb, reflect );
            }
        }
        else
        {
            Vector2 b2 = new Vector2( r + (r-b.x), y ); 

            assert b.x >= a.x;
            if( b.x - a.x > 1.0E-3F )
            {
                Rasterize( bitmap, a, new Vector2( b.x, y ), b.y, rgb, null );
            }
            if( b2.x - a2.x > 1.0E-3F )
            {
                Rasterize( bitmap, a2, b2, b.y, rgb, reflect );
            }
        }
    }

    static private void Rasterize( RawBitmap bitmap, 
                                   Vector2 a,
                                   Vector2 b,
                                   float cy,
                                   int color,
                                   Integer reflect )
    {
        int cr = ((color >> 16) % 256);
        int cg = ((color >> 8) % 256);
        int cb = (color % 256);

        int iflx = (int)Math.round(a.x);
        int ifrx = (int)Math.round(b.x);

        float ml = (b.y - a.y) / (b.x - a.x);
        float mu = (cy - a.y) / (b.x - a.x);

        for( int ix = iflx; ix <= ifrx; ++ix )
        {
            float[] buffer = new float[480];

            float lx = (float)Math.max( ix-0.5F, a.x );
            float rx = (float)Math.min( ix+0.5F, b.x );
            
            float lly = (lx - a.x) * ml + a.y; 
            float lry = (rx - a.x) * ml + a.y;
            
            int illy = (int)Math.round( lly );
            int ilry = (int)Math.round( lry );
            
            float uly = (lx - a.x) * mu + a.y; 
            float ury = (rx - a.x) * mu + a.y;
            
            int iuly = (int)Math.round( uly );
            int iury = (int)Math.round( ury );
               
            for( int y = (int)Math.max( 0, Math.max(illy,ilry) + 1 ); 
                 y <= (int)Math.min( 479, Math.max(iuly,iury) ); ++y )
            {
                buffer[y] = (rx-lx);
            }

            if( ml >= 0 )
            {
                HalfPlane( lx, lly, illy, rx, lry, ilry, ml, ix, 1, buffer );
            }
            else
            {
                HalfPlane( lx, lry, ilry, rx, lly, illy, -ml, ix, 1, buffer );
            }
            if( mu >= 0 )
            {
                HalfPlane( lx, uly, iuly, rx, ury, iury, mu, ix, -1, buffer );
            }
            else
            {
                HalfPlane( lx, ury, iury, rx, uly, iuly, -mu, ix, -1, buffer );
            }

            int tx = reflect != null 
                ? (reflect.intValue() - (ix - reflect.intValue())) : ix;
            for( int y = (int)Math.max( 0, Math.min(illy,ilry) );
                 y <= (int)Math.min( 479, Math.max(iuly,iury) ); ++y )
            {
                float w = buffer[y];
                bitmap.Add( tx, y, (int)(cr*w), (int)(cg*w), (int)(cb*w) );
            }
        }
    }

    static private void HalfPlane( float lx,
                                   float lly,
                                   int illy,
                                   float rx,
                                   float lry,
                                   int ilry,
                                   float ml,
                                   int ix,
                                   int sign,
                                   float[] buffer )
    {
        if( illy == ilry )
        {
            float weight = (rx - lx) * 
                (ilry + 0.5F - 0.5F * lly - 0.5F * lry);
            if( illy >= 0 && illy < 480 )
                buffer[illy] += sign * weight;
        }
        else
        {
            float height = illy + 0.5F - lly;
            float weight1 = 0.5F * height * height / ml;
            if( illy >= 0 && illy < 480 )
                buffer[illy] += sign * weight1;
            
            int y = illy + 1;
            if( illy + 1 < ilry )
            {
                float base = (rx-lx) - height / ml;
                float decr = 1.0F / ml;
                float hdecr = 0.5F*decr;
                base -= 0.5F * decr;
                while( base > hdecr )
                {
                    float weight = (rx-lx) - base;
                    if( y >= 0 && y < 480 )
                        buffer[y] += sign * weight;
                    base -= decr;
                    ++y;
                }
            }
            
            if( y <= ilry )
            {
                float height2 = lry - (ilry-0.5F);
                float weight = (rx-lx) - 0.5F*height2*height2/ml;
                if( ilry >= 0 && ilry < 480 )
                    buffer[ilry] += sign * weight;
            }
        }
    }
}