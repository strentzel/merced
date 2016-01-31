import gameoxide.GO;
import gameoxide.GOImageCanvas;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

class TestConvexPolygon
{
    static public void main( String[] args )
    {
        GO.Initialize( GO.WINDOW640X480 );

/*        Vector2[] points = new Vector2[5];
        points[0] = new Vector2(  1, 0 );
        points[1] = new Vector2(  1, 1 );
        points[2] = new Vector2( -1, 1 );
        points[3] = new Vector2( -1, 0 );
        points[4] = new Vector2(  0,-1 );
*/

        Vector2[] points = new Vector2[7];
        points[0] = new Vector2(  0.003902F,  0.015678F );
        points[1] = new Vector2(  0.217224F,  0.162290F );
        points[2] = new Vector2(  0.782776F,  0.162290F );
        points[3] = new Vector2(  0.996098F,  0.015678F );
        points[4] = new Vector2(  0.996098F, -0.019690F );
        points[5] = new Vector2(  0.500000F, -0.447910F );
        points[6] = new Vector2(  0.003902F, -0.019690F );

        int MULT = 200;
        int XOFF = 320-MULT/2;
        int YOFF = 240;

        for( int i = 0; i < 7; ++i )
        {
            points[i].x = points[i].x * MULT + XOFF;
            points[i].y = points[i].y * MULT + YOFF;;
        }

        ConvexPolygon C = new ConvexPolygon( points, Material.GetDiamond() );
        RawBitmap bitmap = new RawBitmap( 640, 480 );
        for( int d = -10; d <= 10; ++d )
        {
            float angle = (d/4.0F) * (float)Math.PI / 180.0F;
            LightBeam lb = new LightBeam();
            lb.dir = new Vector2( (float)Math.sin(angle), 
                                  -(float)Math.cos(angle) );
            lb.path = new LinkedList<LightBeam.CrossSection>();
            lb.src = new LightBeam.CrossSection();
            //lb.src.start = new Vector2( -0.5F, 0.0F );
            //lb.src.end = new Vector2(    0.5F, 0.0F );
            lb.src.start = new Vector2( 0.217224F*MULT + XOFF, 0.0F + YOFF );
            lb.src.end = new Vector2(   0.500000F*MULT + XOFF, 0.0F + YOFF );
            lb.src.weight = 0.5F;
            lb.weight = 0.5F;
            lb.inside = true;
            List<LightBeam> beams = new LinkedList<LightBeam>();
            beams.add( lb );

            for( int i = 0; i < 10; ++i )
//        while( true )
            {
                System.out.println( "Iteration " + i );
                boolean hasInside = false;
                ListIterator<LightBeam> it = beams.listIterator();
                while( it.hasNext() )
                {
                    if( it.next().inside )
                    {
                        hasInside = true;
                    }
                }

                if( !hasInside ) break;

                it = beams.listIterator();
                List<LightBeam> next = new LinkedList<LightBeam>();
                while( it.hasNext() )
                {
                    LightBeam beam = it.next();
                    if( beam.inside )
                    {
                        List<LightBeam> result = C.Process( beam );
                        ListIterator<LightBeam> it2 = result.listIterator();
                        while( it2.hasNext() )
                        {
                            LightBeam lb2 = it2.next();
                            if( lb2.weight > 0.01 )
                            {
                                next.add( lb2 );
                            }
                        }
                    }
                    else
                    {
                        next.add( beam );
                    }
                }
                beams = next;

                System.out.println();
                it = beams.listIterator();
                while( it.hasNext() )
                {            
                    System.out.println( it.next().toString() );
                }
            }
            System.out.println();
            ListIterator<LightBeam> it0 = beams.listIterator();
            while( it0.hasNext() )
            {            
                System.out.println( it0.next().toString() );
            }

            ListIterator<LightBeam> it = beams.listIterator();
            while( it.hasNext() )
            {            
                it.next().Draw( bitmap );
            }
        }
        GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );
        bitmap.Draw( canvas );
        C.Draw( canvas );
        GO.Draw( 0, 0, canvas );
        GO.Print();
    }
}
