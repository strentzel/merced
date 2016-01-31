import gameoxide.GO;
import gameoxide.GOImageCanvas;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

class XTestRasterizer
{
    static class Shader implements IRasterizer
    {
        public int GetPixel( int x, int y )
        {
            return 0xFFFFFF;
        }
    }

    static ConvexPolygon Generate( Random R, float expLength )
    {
        List<Vector2> points = new LinkedList<Vector2>();
        float baseAngle = R.nextFloat() * 2 * (float)Math.PI; 
        float offsetAngle = 0;
        int expSides = R.nextInt( 8 ) + 3;
        float expAngle = 2 * (float)Math.PI / expSides;
        float expDev = expAngle / 4;
        final float oneDegree = (float)Math.PI / 180.0F;

        points.add( new Vector2( R.nextInt(320) + 160,
                                 R.nextInt(240) + 120 ) );
        boolean init = true;
        while( true )
        {
            offsetAngle += Math.max( 
                (float)R.nextGaussian() * expDev + expAngle, oneDegree );

            Vector2 last = points.get( points.size()-1 );
            
            float length = Math.max( 
                (float)R.nextGaussian() * expLength/4.0F + expLength, 
                1.0F );
            float angle = baseAngle + offsetAngle;
            Vector2 next = new Vector2( 
                last.x + length * (float)Math.cos( angle ),
                last.y + length * (float)Math.sin( angle ) );
            
            if( !init )
            {
                Vector2 first = points.get( 0 );
                Vector2 second = points.get( 1 );
                Vector2 edge = new Vector2();
                Vector2.Minus( edge, next, last );
                Vector2 close = new Vector2();
                Vector2.Minus( close, first, last );
                Vector2 future = new Vector2();
                Vector2.Minus( future, first, next );
                Vector2 start = new Vector2();
                Vector2.Minus( start, second, first );

                if( edge.Cross( close ) <= 0 ||
                    future.Cross( start ) <= 0 )
                {
                    break;
                }
            }
            else
            {
                init = false;
            }

            points.add( next );
        }

        Vector2[] array = new Vector2[points.size()];
        ListIterator<Vector2> iter = points.listIterator();
        int i = 0;
        while( iter.hasNext() )
        {
            array[i++] = iter.next();
        }
        return new ConvexPolygon( array );
    }

    static public void main( String[] args )
    {
        GO.Initialize( GO.WINDOW640X480 );
        GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );

        Random R = new Random( 117 );
        Shader S = new Shader();

        while(true)
        {
            ConvexPolygon p = Generate( R, 100 );
            RawBitmap bitmap = new RawBitmap( 640, 480 );
            p.Rasterize( bitmap, 0xFF0000 );

            bitmap.Draw( canvas );
            GO.Draw( 0, 0, canvas );
            GO.Print();

            try
            {
                Thread.sleep( 1000 );
            }
            catch( Exception ex )
            {
            }

            p.Rasterize( bitmap, S );
            bitmap.Draw( canvas );
            GO.Draw( 0, 0, canvas );
            GO.Print();

            try
            {
                Thread.sleep( 4000 );
            }
            catch( Exception ex )
            {
            }
        }
    }
 }
