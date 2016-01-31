import gameoxide.GO;
import gameoxide.GOImageCanvas;
import java.util.Random;

class XTestConvexPolyhedra
{
    static public void main( String[] args )
    {
        GO.Initialize( GO.WINDOW640X480 );
        GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );

        ConvexPolyhedra P = ConvexPolyhedra.MakeTestCube();
        Random R = new Random();

        while(true)
        {
            P.Transform(Matrix3.MakeXRotation(R.nextFloat()*2*(float)Math.PI));
            P.Transform(Matrix3.MakeYRotation(R.nextFloat()*2*(float)Math.PI));
            P.Transform(Matrix3.MakeZRotation(R.nextFloat()*2*(float)Math.PI));

            System.out.println( "GO" );

            RawBitmap bitmap = new RawBitmap( 640, 480 );
            P.Rasterize( 100.0F, new Vector2( 320, 240 ), bitmap, 0x00FFFF );
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
        }
    }
}
