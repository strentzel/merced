import gameoxide.GO;
import gameoxide.GOImageCanvas;

class XTestDiamond
{
    static public void main( String[] args )
    {
        GO.Initialize( GO.WINDOW640X480 );
        GOImageCanvas canvas = GO.GetImageCanvas( 640, 480 );

        RawBitmap bitmap = new RawBitmap( 640, 480 );

        RoundCut diamond = RoundCut.MakeTest();
        System.out.println( "VLM=" + diamond.ComputeVolume() );
        Matrix3 rot = new Matrix3();
        Matrix3.Multiply( rot, 
                          Matrix3.MakeYRotation( (float)(Math.PI/180.0*60) ),
                          Matrix3.MakeXRotation( (float)(Math.PI/180.0*45) ) );

        Matrix3 irot = new Matrix3();
        irot.d00 =  400.0F;
        irot.d10 =    0.0F;
        irot.d20 =    0.0F;
        irot.d01 =    0.0F;
        irot.d11 =  400.0F;
        irot.d21 =    0.0F;
        irot.d02 =    0.0F;
        irot.d12 =    0.0F;
        irot.d22 =  400.0F;

        diamond.Transform( diamond.GetCrossSectionRot(0) );
        diamond.Rasterize( 400.0F, new Vector2( 320, 240 ), bitmap, 0x007F7F );

        Vector3[] cs = diamond.GetCrossSection(0);
        for( int i = 0; i < cs.length-1; ++i )
        {
            int in = (i+1)%cs.length;
            Vector2 a = new Vector2( 400*cs[i].x+320, 400*cs[i].y+240 );
            Vector2 b = new Vector2( 400*cs[in].x+320, 400*cs[in].y+240 );
            Line l = new Line( a, b, 5 );
            l.Rasterize( bitmap, 0xFF0000 );
        }

        bitmap.Draw( canvas );
        GO.Draw( 0, 0, canvas );
        GO.Print();

        try
        {
            Thread.sleep( 1000 );
        }
        catch( Exception ex ){}
        GO.Print();
    }
}
