import gameoxide.GO;
import gameoxide.GOImageCanvas;

class TestTriangle
{
    public static void main( String[] args )
    {

        if( args.length == 2 )
        {
            int h0 = Integer.parseInt( args[0] );
            int h1 = Integer.parseInt( args[1] );

            RawBitmap bitmap = new RawBitmap( 640, 480 );
            Triangle t = new Triangle( new Vector2( 100, 50 ),
                                       new Vector2( 105, 50+h0 ),
                                       new Vector2( 110, 50+h1 ) );
            float theta = 0.1F;
            float r = 100;
            float incr = (float)(Math.PI*2/3.0);
            t = new Triangle( 
                new Vector2( (float)(320 + r*Math.cos(theta)),
                             (float)(240 + r*Math.sin(theta)) ),
                new Vector2( (float)(320 + r*Math.cos(theta+incr)),
                             (float)(240 + r*Math.sin(theta+incr)) ),
                new Vector2( (float)(320 + r*Math.cos(theta+2*incr)),
                             (float)(240 + r*Math.sin(theta+2*incr)) ) );
            t = new Triangle(
                new Vector2( 320, 250 ),
                new Vector2( 325, 240 ),
                new Vector2( 320, 260 ) );

            if( h0 == 1 )
            {
                t = new Triangle( 
                    new Vector2( 301.13266F, 240.0F ),
                    new Vector2( 320.0F, 150.41798F ),
                    new Vector2( 320.0F, 240.0F ) );
            }
            else
            {
                t = new Triangle( 
                    new Vector2( 301.13266F, 166.70418F ),
                    new Vector2( 301.13266F, 240.0F ),
                    new Vector2( 320.0F, 150.41798F ) );
            }
            t = new Triangle(
                new Vector2( 263.44482F, 240.0F ),
                new Vector2( 267.3026F, 195.90514F ),
                new Vector2( 320.0F, 150.418F ) );
            t.Rasterize( bitmap, 0xBFBFBF );

            int mxx = 320;
            int mnx = 320;
            int mxy = 240;
            int mny = 240;
            for( int x = 0; x < 640; ++x )
            {
                for( int y = 0; y < 480; ++y )
                {
                    if( bitmap.data[y*640+x] > 0 )
                    {
                        if( x > mxx ) mxx = x;
                        if( x < mnx ) mnx = x;
                        if( y > mxy ) mxy = y;
                        if( y < mny ) mny = y;
                    }
                }
            }

            for( int y = mxy; y >= mny; --y )
            {
                for( int x = mnx; x <= mxx; ++x )
                {
                    System.out.format( "%5d ", (int)(bitmap.data[y*640+x])>>16 );
                }
                System.out.println();
            }
        }
        else
        {
            GO.Initialize( GO.WINDOW640X480 );
            GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );

            float incr = (float) (2 * Math.PI / 3.0F);
            float theta = 0;
            float r = 100;
            float rate = (float) (2 * Math.PI / (30 * 60));
            long sum = 0;
            int cnt = 0;
            while( true )
            {
                RawBitmap bitmap = new RawBitmap( 640, 480 );
                long start = System.currentTimeMillis();
                Triangle t = new Triangle( 
                    new Vector2( (float)(320 + r*Math.cos(theta)),
                                 (float)(240 + r*Math.sin(theta)) ),
                    new Vector2( (float)(320 + r*Math.cos(theta+incr)),
                                 (float)(240 + r*Math.sin(theta+incr)) ),
                    new Vector2( (float)(320 + r*Math.cos(theta+2*incr)),
                                 (float)(240 + r*Math.sin(theta+2*incr)) ) );
//                t.Rasterize( bitmap, 0xFFFFFF );
                long stop = System.currentTimeMillis();
                sum += (stop-start);
                cnt++;
                System.out.println( sum + "/" + cnt + "=" + (sum/cnt) ); 
/*
                for( int i = 0; i < 640*480; ++i )
                {
                    int w = bitmap.data[i] * 255 / 10000;
                    if( w > 255 ) w = 255;
                    bitmap.data[i] = (w << 16) + (w << 8) + (w);
                }
*/
                Triangle t1 = new Triangle( 
                    new Vector2( 301.13266F, 240.0F ),
                    new Vector2( 320.0F, 150.41798F ),
                    new Vector2( 320.0F, 240.0F ) );
                Triangle t2 = new Triangle( 
                    new Vector2( 301.13266F, 166.70418F ),
                    new Vector2( 301.13266F, 240.0F ),
                    new Vector2( 320.0F, 150.41798F ) );

                t1.Rasterize( bitmap, 0xBFBFBF );
                t2.Rasterize( bitmap, 0xBFBFBF );

                bitmap.Draw( canvas );
                GO.Draw( 0, 0, canvas );
                GO.Print();

                theta += rate;
                while( theta > Math.PI * 2 )
                    theta -= Math.PI * 2;
            }
        }            
    }
}
