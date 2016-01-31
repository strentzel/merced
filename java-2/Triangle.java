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
        System.out.println( "TRI: " + a + " " + b + " " + c + " " + (rgb>>16) );
        float y = (b.x-a.x)*(c.y-a.y)/(c.x-a.x) + a.y;
        float r = (float)Math.round( c.x );
        Vector2 a2 = new Vector2( r + (r-c.x), c.y );
        if( b.y < y )
        {
            Vector2 b2 = new Vector2( r + (r-b.x), b.y ); 

            assert b.x >= a.x;
            if( b.x - a.x > 1.0E-3F )
            {
                Rasterize( bitmap, a, b, y, rgb, false );
            }
            if( b2.x - a2.x > 1.0E-3F )
            {
                Rasterize( bitmap, a2, b2, y, rgb, true );
            }
        }
        else
        {
            Vector2 b2 = new Vector2( r + (r-b.x), y ); 

            assert b.x >= a.x;
            if( b.x - a.x > 1.0E-3F )
            {
                Rasterize( bitmap, a, new Vector2( b.x, y ), b.y, rgb, false );
            }
            if( b2.x - a2.x > 1.0E-3F )
            {
                Rasterize( bitmap, a2, b2, b.y, rgb, true );
            }
        }
    }


    static public void Rasterize( RawBitmap bitmap, 
                                  Vector2 a,
                                  Vector2 b,
                                  float cy,
                                  int color,
                                  boolean reflect )
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
               
            for( int y = (int)Math.max(illy,ilry) + 1; 
                 y <= (int)Math.max(iuly,iury); ++y )
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

            int tx = reflect ? (iflx - (ix - iflx)) : ix;
            for( int y = (int)Math.min(illy,ilry);
                 y <= (int)Math.max(iuly,iury); ++y )
            {
                float w = buffer[y];
                bitmap.Add( tx, y, (int)(cr*w), (int)(cg*w), (int)(cb*w) );
            }
        }
    }

    static public void HalfPlane( float lx,
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
            buffer[illy] += sign * weight;
        }
        else
        {
            float height = illy + 0.5F - lly;
            float weight1 = 0.5F * height * height / ml;
            buffer[illy] += sign * weight1;
            
            int y = illy + 1;
            if( illy + 1 < ilry )
            {
                float base = (rx-lx) - height / ml;
                float decr = 1.0F / ml;
                base -= 0.5F * decr;
                while( base > 0 )
                {
                    float weight = (rx-lx) - base;
                    buffer[y] += sign * weight;
                    base -= decr;
                    ++y;
                }
            }
            
            if( y <= ilry )
            {
                float height2 = lry - (ilry-0.5F);
                float weight = (rx-lx) - 0.5F*height2*height2/ml;
                buffer[ilry] += sign * weight;
            }
        }
    }
/*
    static public void Rasterize( RawBitmap bitmap, 
                                   Vector2 a,
                                   Vector2 b,
                                   float cy )
    {
        int iflx = (int)Math.round(a.x);
        int ifrx = (int)Math.round(b.x);

        float ml = (b.y - a.y) / (b.x - a.x);
        float mu = (cy - a.y) / (b.x - a.x);

        for( int ix = iflx; ix <= ifrx; ++ix )
        {
            float lx = (float)Math.max( ix-0.5F, a.x );
            float rx = (float)Math.min( ix+0.5F, b.x );

            float lly = (lx - a.x) * ml + a.y; 
            float lry = (rx - a.x) * ml + a.y;

            int illy = (int)Math.round( lly );
            int ilry = (int)Math.round( lry );
            
            if( illy == ilry )
            {
                float weight = (rx - lx) * 
                    (ilry + 0.5F - 0.5F * lly - 0.5F * lry);
                bitmap.data[illy * bitmap.width + ix] = (int)(10000*weight);
            }
            else
            {
                float height = illy + 0.5F - lly;
                float weight1 = 0.5F * height * height / ml;
                bitmap.data[illy * bitmap.width + ix] = (int)(10000*weight1);

                if( illy + 1 < ilry )
                {
                    float base = (rx-lx) - height / ml;
                    float decr = 1.0F / ml;
                    base -= 0.5F * decr;
                    int y = illy + 1;
                    while( base > 0 )
                    {
                        float weight = (rx-lx) - base;
                        bitmap.data[y * bitmap.width + ix] = (int)(10000*weight);
                        base -= decr;
                        ++y;
                    }
                }
                
                float height2 = lry - (ilry-0.5F);
                float weight = (rx-lx) - 0.5F*height2*height2/ml;
                bitmap.data[ilry * bitmap.width + ix] = (int)(10000*weight);
            }
        }
    }
*/

/*
    static public void RasterizeLeft( RawBitmap bitmap, 
                                      Vector2 a,
                                      Vector2 b,
                                      float cy )
    {
        int iflx = (int)Math.round(a.x);
        int ifrx = (int)Math.round(b.x);

        float ml = (b.y - a.y) / (b.x - a.x);
        float mu = (cy - a.y) / (b.x - a.x);

        for( int ix = iflx; ix <= ifrx; ++ix )
        {
            int[] buffer = new int[480];

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
                       
            //for( int y = ilry + 1; y <= iury; ++y )
            for( int y = (int)Math.max(illy,ilry) + 1; 
                 y <= (int)Math.max(iuly,iury); ++y )
            {
                buffer[y] = (int)((rx-lx)*10000);
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

//            for( int y = illy; y <= iury; ++y )
            for( int y = (int)Math.min(illy,ilry);
                 y <= (int)Math.max(iuly,iury); ++y )
            {
                bitmap.data[y*bitmap.width + ix] += buffer[y];
            }
        }
    }

    static public void RasterizeRight( RawBitmap bitmap, 
                                       Vector2 a,
                                       Vector2 b,
                                       float cy )
    {
        int iflx = (int)Math.round(a.x);
        int ifrx = (int)Math.round(b.x);

        int ceil = (int)Math.ceil(a.x);
        int floor = (int)Math.floor(a.x);

        float ml = (b.y - a.y) / (b.x - a.x);
        float mu = (cy - a.y) / (b.x - a.x);

        for( int ix = iflx; ix <= ifrx; ++ix )
        {
            int[] buffer = new int[480];

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
                       
            //for( int y = ilry + 1; y <= iury; ++y )
            for( int y = (int)Math.max(illy,ilry) + 1; 
                 y <= (int)Math.max(iuly,iury); ++y )
            {
                buffer[y] = (int)((rx-lx)*10000);
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

            int tx = iflx - (ix - iflx);
//            for( int y = illy; y <= iury; ++y )
            for( int y = (int)Math.min(illy,ilry);
                 y <= (int)Math.max(iuly,iury); ++y )
            {
                bitmap.data[y*bitmap.width + tx] += buffer[y];
            }
        }
    }

    static public void HalfPlaneD( float lx,
                                   float lly,
                                   int illy,
                                   float rx,
                                   float lry,
                                   int ilry,
                                   float ml,
                                   int ix,
                                   int sign,
                                   int[] buffer )
    {
        if( illy == ilry )
        {
            float weight = (rx - lx) * 
                (ilry + 0.5F - 0.5F * lly - 0.5F * lry);
            buffer[illy] += (int)(10000*weight);
        }
        else
        {
            float height = ilry + 0.5F - lry;
            float weight1 = 0.5F * height * height / -ml;
            buffer[ilry] += sign * (int)(10000*weight1);

            if( ilry + 1 < illy )
            {
                float base = (rx-lx) - height / -ml;
                float decr = 1.0F / -ml;
                base -= 0.5F * decr;
                int y = ilry + 1;
                while( base > 0 )
                {
                    float weight = (rx-lx) - base;
                    buffer[y] += sign * (int)(10000*weight);
                    base -= decr;
                    ++y;
                }
            }
            
            float height2 = lly - (illy-0.5F);
            float weight = (rx-lx) - 0.5F*height2*height2/-ml;
            buffer[illy] += sign * (int)(10000*weight);
        }
    }
*/
}