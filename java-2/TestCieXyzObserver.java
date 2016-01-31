class TestCieXyzObserver
{
    static public void main( String[] args )
    {
        CieXyzObserver observer = null;
        try
        {
            observer = new CieXyzObserver( "data/cie.csv" );
        }
        catch( java.io.IOException ex ){System.exit(0);}

        for( int i = 0; i <= 100; ++i )
        {
            float w = i/100.0F;
            CieXyz xyz = new CieXyz();
            xyz = xyz.AddMult( observer.Integrate( 0, 1000 ), w );
            int rgb = xyz.toRgb( 0x000000 );
            float r = ((rgb>>16)%256)/255.0F;
            float g = ((rgb>>8)%256)/255.0F;
            float b = (rgb%256)/255.0F;
            System.out.println( w + " " + r + " " + g + " " + b );
        }
    }
}