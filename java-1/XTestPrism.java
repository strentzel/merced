import gameoxide.GO;
import gameoxide.GOImageCanvas;
import gameoxide.GOInput;
import java.util.Random;

class XTestPrism
{
    static public void main( String[] args )
        throws java.lang.Exception
    {
        GO.Initialize( GO.WINDOW640X480 );
        GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );

        Random R = new Random();

        Material material = Material.GetDiamond();
        CieXyzObserver observer = new CieXyzObserver( "cie.csv" );

        Vector2[] points = new Vector2[4];
        points[0] = new Vector2( 280, 210 ); 
        points[1] = new Vector2( 360, 210 ); 
        points[2] = new Vector2( 360, 270 ); 
        points[3] = new Vector2( 280, 270 ); 
        ConvexPolygon centerBox = new ConvexPolygon( points );

        while( true )
        {
            Vector2 surface1 = new Vector2( 640 * R.nextFloat(),
                                            480 * R.nextFloat() );
            float angle = 2 * (float)Math.PI * R.nextFloat();
            float length = (float) Math.min( 
                250, 50 * Math.exp( R.nextGaussian() ) );
            Vector2 surface2 = new Vector2( 
                surface1.x + (float)Math.cos( angle ) * length,
                surface1.y + (float)Math.sin( angle ) * length );
            Vector2 center = new Vector2(
                surface1.x + (float)Math.cos( angle ) * length / 2,
                surface1.y + (float)Math.sin( angle ) * length / 2 );

            for( int i = 0; i < 25; ++i )
            {
                Vector2 mouse = new Vector2( GOInput.GetX(), GOInput.GetY() );
                Vector2 lightDir = new Vector2();
                Vector2.Minus( lightDir, center, mouse );
                Vector2.Normalize( lightDir );
                
                Prism P = new Prism( material, observer, surface1, surface2,
                                     lightDir, 1.0F );

                canvas.setColor( 0x000000 );
                canvas.fill( 0, 0, 640, 480 );

                RawBitmap bitmap = new RawBitmap( 640, 480 );
                ConvexPolygon intersect = P.Intersect( centerBox );
                centerBox.Rasterize( bitmap, 0xFF0000 );
                if( intersect != null )
                    intersect.Rasterize( bitmap, 0x00FF00 );

                bitmap.Draw( canvas );

                Vector2 sp2 = new Vector2();
                Vector2.Add( sp2, P.mSurfacePoint, P.mSurfaceVec );
                Vector2 mn = new Vector2();
                Vector2.Scale( mn, P.mMinDir, 50 );
                Vector2 mn1 = new Vector2();
                Vector2.Add( mn1, P.mSurfacePoint, mn );
                Vector2 mn2 = new Vector2();
                Vector2.Add( mn2, sp2, mn );
                Vector2 mx = new Vector2();
                Vector2.Scale( mx, P.mMaxDir, 50 );
                Vector2 mx1 = new Vector2();
                Vector2.Add( mx1, P.mSurfacePoint, mx );
                Vector2 mx2 = new Vector2();
                Vector2.Add( mx2, sp2, mx );
                canvas.setColor( 0xFFFFFF );
                canvas.moveTo( (int)P.mSurfacePoint.x,
                               (int)P.mSurfacePoint.y );
                canvas.drawTo( (int)sp2.x, (int)sp2.y );
                canvas.setColor( 0x00FFFF );
                canvas.moveTo( (int)P.mSurfacePoint.x,
                               (int)P.mSurfacePoint.y );
                canvas.drawTo( (int)mn1.x, (int)mn1.y );
                canvas.moveTo( (int)sp2.x,(int)sp2.y );
                canvas.drawTo( (int)mn2.x, (int)mn2.y );
                canvas.setColor( 0xFF0000 );
                canvas.moveTo( (int)P.mSurfacePoint.x,
                               (int)P.mSurfacePoint.y );
                canvas.drawTo( (int)mx1.x, (int)mx1.y );
                canvas.moveTo( (int)sp2.x,(int)sp2.y );
                canvas.drawTo( (int)mx2.x, (int)mx2.y );

                GO.Draw( 0, 0, canvas );
                GO.Print();

                Thread.sleep( 200 );
            }
        }
    }
}
