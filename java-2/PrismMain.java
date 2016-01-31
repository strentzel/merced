import gameoxide.*;

public class PrismMain {

	static final float ThetaD = 20.0F*(float)Math.PI/180;
	
	static float GetWeight( int rgb )
	{
	    int r = (rgb >> 16) & 0xFF;
	    int g = (rgb >> 8) & 0xFF;
	    int b = rgb & 0xFF;
	    return (r+g+b)/3.0F;
	}
	
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
		
		float minN = 4834000.0F / 2008000.0F;
		float minTheta = (float)Math.asin(minN*Math.sin(ThetaD));
		float maxN = (float)(((731.3-380)*(731.3-380)+4834000.0)/2008000.0);
		float maxTheta = (float)Math.asin(maxN*Math.sin(ThetaD));
		
		int light = 12;
		
		CieXyz best = new CieXyz();
		float maxY = 0;
		float maxW = 0;
		
		GOImageCanvas canvas = GO.GetImageCanvas(640,480);
		for( int row = 0; row < 480; ++row )
		{
			System.out.println( "Row " + row + " " + maxW );
			for( int col = 0; col < 640; ++col )
			{
				CieXyz xyz = new CieXyz();

				float thetaLow;
				float thetaHigh;
				int ls;
				int rs;
				{
				    float x = 320 - light / 2.0F + 0 / 10.0F + 0.05F;
				    float distance = (float)Math.sqrt((col-x)*(col-x)+(480-row)*(480-row));
				    float theta = (float)Math.atan(Math.abs(col-x)/(float)(480-row));
				    float theta0 = theta - (float)Math.asin(0.5/distance);
				    float theta1 = theta + (float)Math.asin(0.5/distance);
				    thetaLow = theta0;
				    thetaHigh = theta1;
				    if( col >= x ) ls = 1; else ls = -1;
				}
				{
				    float x = 320 - light / 2.0F + (10*light-1) / 10.0F + 0.05F;
				    float distance = (float)Math.sqrt((col-x)*(col-x)+(480-row)*(480-row));
				    float theta = (float)Math.atan(Math.abs(col-x)/(float)(480-row));
				    float theta0 = theta - (float)Math.asin(0.5/distance);
				    float theta1 = theta + (float)Math.asin(0.5/distance);
				    thetaLow = (float)Math.min( thetaLow, theta0 );
				    thetaHigh = (float)Math.max( thetaHigh, theta1 );
				    if( col > x ) rs = 1; else rs = -1;
				}
				
				if( thetaLow < maxTheta && thetaHigh > minTheta ||
				    ls == 1 && rs == -1 )
				{
                    for( int i = 0; i < 10 * light; ++i )
                    {
                        float x = 320 - light / 2.0F + i / 10.0F + 0.05F;
                        float distance = (float)Math.sqrt((col-x)*(col-x)+(480-row)*(480-row));
                        float theta = (float)Math.atan(Math.abs(col-x)/(float)(480-row));
                        float theta0 = theta - (float)Math.asin(0.5/distance);
                        float theta1 = theta + (float)Math.asin(0.5/distance);
                        if( (theta0 < maxTheta || theta1 > minTheta) && col >= x )
                        {
                            float n0 = (float)Math.max(minN,(Math.sin(theta0)/Math.sin(ThetaD)));
                            float n1 = (float)Math.min(maxN,(Math.sin(theta1)/Math.sin(ThetaD)));
        
                            float w0 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n0 - 4834000,0)));
                            float w1 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n1 - 4834000,0)));
                            
                            xyz = xyz.AddMult( observer.Integrate( w1, w0 ), 0.1F );
                        }
                    }
                }
                else
                {
                    //System.out.println( "Tl " + thetaLow + " Th " + thetaHigh + " Tn " + minTheta + " Tx " + maxTheta );
                }
                int rgb = xyz.toRgb( 0x0 );
                
                float w = GetWeight( rgb );
                if( w > 250 ){
                    rgb = 0xFF0000;
                }
                else if( w > 225 ) rgb = 0x00FF00;
                else if( w > 200 )
                {
                    //System.out.println( "rgb " + (rgb>>16) + " " + ((rgb>>8)&0xff) + " " + (rgb&0xff) + "Tl " + thetaLow + " Th " + thetaHigh + " Tn " + minTheta + " Tx " + maxTheta );
                    rgb = 0x0000FF;
                }
                else if( w > 0  ) rgb = 0x00FFFF;
                rgb = xyz.toRgb( 0x0 );
                
                canvas.setColor( rgb );
                
				if( xyz.X > maxY )
				{
					maxY = xyz.X;
					best = xyz;
				}
				if( xyz.Y > maxY )
				{	
					maxY = xyz.Y;
					best = xyz;
				}
				if( xyz.Z > maxY )
				{
					maxY = xyz.Z;
					best = xyz;
				}
				maxW = (float)Math.max(maxW,GetWeight( xyz.toRgb(0x0) ) );

/*
					float k = (w0-w1)/100.0F;
					for( int i = 0; i < 100; ++i )
					{
						float n = n0 + (n1-n0)*i/100.0F;
						float wavelength = (float)(731.3 - Math.sqrt(Math.max(2008000 * n - 4834000,0)));
						datum = datum.Add( observer.GetDatum( wavelength ).Multiply( k ) );
					}
					datum = datum.Multiply( 1.0F / 21.3715F / 5.0F );
					float s = 1.0F;//(float)((n1-n0)/(maxN-minN));
					xyz.X = datum.x * s; if( xyz.X > 1 ) xyz.X = 1;
					xyz.Y = datum.y * s; if( xyz.Y > 1 ) xyz.Y = 1;
					xyz.Z = datum.z * s; if( xyz.Z > 1 ) xyz.Z = 1;
					canvas.setColor( xyz.toRgb() );
					
					if( datum.z > maxY )
					{
						maxY = datum.z;
						best = datum;
					}
*/

					//canvas.setColor( (int)(255*(theta+Math.PI)/Math.PI) );
/*				
				if( 0 == xyz.toRgb( 0x0 ) || 0xffffff == xyz.toRgb( 0x0 ) )
				{
				    canvas.setColor( 0x0 );
				}
				else
				{
				    canvas.setColor( 0x00ffff );
				}
				*/
				canvas.moveTo(col, row);
				canvas.drawTo(col, row+1 );
				
				//float s = row/480.0F;//distance / 576.0F;
				//float wavelength = (380 + (639-col)/640.0F * (780-380));
				//float s = (wavelength-380)/400/2 + 0.5F;
				//(float)java.lang.Math.exp((row/480.0F)*20-10);
				//float s = row / 480.0F; 

//				float x = col / 640.0F;
//				float y = row / 480.0F;
//				xyz.Y = 0.25F;
//				xyz.X = xyz.Y*x/y;
//				xyz.Z = xyz.Y*(1-x-y)/y;
			}
		}
	
		System.out.println("Ready");
		GO.Draw(0, 0, canvas);
		GO.Print();
	}

}
