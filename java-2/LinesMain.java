import gameoxide.*;

public class LinesMain {

    /**
     * @param args
     */
    public static void main(String[] args) {
        CieXyzObserver observer;
	
        try
        {
            observer = new CieXyzObserver( "data/cie.csv" );
        }
        catch( java.io.IOException e )
        {
            System.out.println( e.toString() );
            return;
        }
	
        // TODO Auto-generated method stub
        GO.Initialize(GO.WINDOW640X480);
        GOImageCanvas canvas = GO.GetImageCanvas(640,480);
        
        Vector2 tmp1 = new Vector2();
        Vector2 tmp2 = new Vector2();
        Vector2 tmp3 = new Vector2();

        while( true )
        {
            Prism P = new Prism( observer,
                                 new Vector2( 320, 240 ),
//                                 new Vector2( 450, 263 ),
//                                 new Vector2( 429, 287 ),
                                 new Vector2( GOInput.GetX(),
                                              GOInput.GetY() ),
                                 new Vector2( 0, 1 ) );
            
            for( int row = 0; row < 480; ++row )
            {
                tmp1.y = row;
                for( int col = 0; col < 640; ++col )
                {
                    tmp1.x = col;
                    canvas.setColor( P.Compute( tmp1 ) );
                    canvas.moveTo( col, row );
                    canvas.drawTo( col, row+1 );
                }
            }
//            canvas.fill( 0, 0, 640, 480 );

            

            canvas.setColor( 0xFFFFFF );
            canvas.moveTo((int)P.mSurfacePoint.x, (int)P.mSurfacePoint.y);
            Vector2.Add( tmp1, P.mSurfacePoint, P.mSurfaceVec );
            canvas.drawTo((int)tmp1.x, (int)tmp1.y);

            canvas.setColor( 0x00FFFF );
            canvas.moveTo((int)P.mSurfacePoint.x, (int)P.mSurfacePoint.y);
            Vector2.Scale( tmp2, P.mMinDir, 100 );
            Vector2.Add( tmp3, P.mSurfacePoint, tmp2 );
            canvas.drawTo((int)tmp3.x, (int)tmp3.y);

            canvas.setColor( 0xFF0000 );
            canvas.moveTo((int)P.mSurfacePoint.x, (int)P.mSurfacePoint.y);
            Vector2.Scale( tmp2, P.mMaxDir, 100 );
            Vector2.Add( tmp3, P.mSurfacePoint, tmp2 );
            canvas.drawTo((int)tmp3.x, (int)tmp3.y);

            canvas.setColor( 0x00FFFF );
            canvas.moveTo((int)tmp1.x, (int)tmp1.y);
            Vector2.Scale( tmp2, P.mMinDir, 100 );
            Vector2.Add( tmp3, tmp1, tmp2 );
            canvas.drawTo((int)tmp3.x, (int)tmp3.y);

            canvas.setColor( 0xFF0000 );
            canvas.moveTo((int)tmp1.x, (int)tmp1.y);
            Vector2.Scale( tmp2, P.mMaxDir, 100 );
            Vector2.Add( tmp3, tmp1, tmp2 );
            canvas.drawTo((int)tmp3.x, (int)tmp3.y);

//            System.out.println( "X " + GOInput.GetX() + " Y " + GOInput.GetY() );

            GO.Draw(0, 0, canvas);
            GO.Print();
	}
    }
}
