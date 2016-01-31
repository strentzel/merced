public class Prism {

    static public final float MinN = 4834000.0F / 2008000.0F;
    static public final float MaxN = (float)(((731.3-380)*(731.3-380)+4834000.0)/2008000.0);
    static public final float Ninety = (float)(Math.PI / 2.0);
    static public final double ExSinTheta = Math.cos( 0.0001 );

    CieXyzObserver mObserver;
    public Vector2 mSurfacePoint;
    public Vector2 mSurfaceVec;
    public Vector2 mLightDir;
    public float mThetaD;
    public float mMinTheta;
    public float mMaxTheta;
    public Vector2 mMinDir;
    public Vector2 mMaxDir;
    float mOrigWeight;

    Vector2 mNegB;
    Vector2 mNegC;

    Vector2 mNormA;
    Vector2 mNormB;
    Vector2 mNormC;
    Vector2 mNormD;
    Vector2 mNormE;

    int mCoreRgb;

    Vector2 mTmp1;
    Vector2 mTmp2;

    public Prism( CieXyzObserver obs,
                  Vector2 surface1,
                  Vector2 surface2,
                  Vector2 lightDir,
                  float origWeight )
    {
        mOrigWeight = origWeight;
        mTmp1 = new Vector2();
        mTmp2 = new Vector2();

        Vector2 tmp1 = new Vector2();
        Vector2 tmp2 = new Vector2();

        mObserver = obs;
        mSurfaceVec = new Vector2();
        Vector2 surfaceDir;
        mLightDir = new Vector2( lightDir );

        Vector2.Minus( tmp1, surface2, surface1 );
        if( mLightDir.Dot( tmp1 ) >= 0.0F )
        {
            mSurfacePoint = new Vector2( surface1 );
            Vector2.Minus( mSurfaceVec, surface2, surface1 );
        }
        else
        {
            mSurfacePoint = new Vector2( surface2 );
            Vector2.Minus( mSurfaceVec, surface1, surface2 );
        }
        Vector2.Minus( tmp1, surface2, surface1 );
        if( tmp1.Norm2() > 1e-8F )
        {
            surfaceDir = new Vector2( mSurfaceVec );
            Vector2.Normalize( surfaceDir );
        }
        else
        {
            mSurfacePoint = new Vector2( surface1 );
            surfaceDir = new Vector2( mLightDir );
            Vector2.Normalize( surfaceDir );
        }
        mThetaD = Ninety - 
            (float) Math.acos( surfaceDir.Dot( mLightDir ) /
                               mLightDir.Norm());
        Vector2 normal = new Vector2( surfaceDir.y, -surfaceDir.x );
        Vector2.Normalize( normal );
        if( mLightDir.Dot( normal ) < 0 )
        {
            Vector2.Scale( normal, normal, -1 );
        }

        mNegB = new Vector2();
        mNegC = new Vector2();
        mNormA = new Vector2();
        mNormB = new Vector2();
        mNormC = new Vector2();
        mNormD = new Vector2();
        mNormE = new Vector2();

        double minSinTheta = MinN*Math.sin(mThetaD);
        if( -ExSinTheta < minSinTheta && minSinTheta < ExSinTheta )
        {
            mMinTheta = (float)Math.asin(MinN*Math.sin(mThetaD));
            mMinDir = new Vector2();
            Vector2.Scale( tmp1, surfaceDir, (float)Math.sin( mMinTheta ) );
            Vector2.Scale( tmp2, normal, (float)Math.cos( mMinTheta ) );
            Vector2.Add( mMinDir, tmp1, tmp2 );
            Vector2.Normalize( mMinDir );

            Vector2.Scale( tmp1, surfaceDir, -1 );
            Vector2.NormFromRel( mNormA, mMinDir, surfaceDir );
            Vector2.NormFromRel( mNormC, mMinDir, tmp1 );
            Vector2.Scale( mNegC, mMinDir, -1/mMinDir.Norm() );
        }
        else
        {
            mMinTheta = Ninety;
            mMinDir = new Vector2( surfaceDir );
            Vector2.Scale( tmp1, lightDir, -1 );
            Vector2.NormFromRel( mNormA, mMinDir, tmp1 );
            Vector2.NormFromRel( mNormC, mMinDir, lightDir );
        }

        double maxSinTheta = MaxN*Math.sin(mThetaD);
        if( -ExSinTheta < maxSinTheta && maxSinTheta < ExSinTheta )
        {
            mMaxTheta = (float)Math.asin(MaxN*Math.sin(mThetaD));
            mMaxDir = new Vector2();
            Vector2.Scale( tmp1, surfaceDir, (float)Math.sin( mMaxTheta ) );
            Vector2.Scale( tmp2, normal, (float)Math.cos( mMaxTheta ) );
            Vector2.Add( mMaxDir, tmp1, tmp2 );
            Vector2.Normalize( mMaxDir );
 
            Vector2.Scale( tmp1, surfaceDir, -1 );
            Vector2.NormFromRel( mNormB, mMaxDir, surfaceDir );
            Vector2.NormFromRel( mNormD, mMaxDir, tmp1 );
            Vector2.Scale( mNegB, mMaxDir, -1/mMaxDir.Norm() );
        }
        else
        {
            mMaxTheta = Ninety;
            mMaxDir = new Vector2( surfaceDir );

            Vector2.Scale( tmp1, lightDir, -1 );
            Vector2.NormFromRel( mNormB, mMaxDir, tmp1 );
            Vector2.NormFromRel( mNormD, mMaxDir, lightDir );
        }

        if( surfaceDir.Dot( mLightDir ) > 1E-8 )
        {
            Vector2.NormFromRel( mNormE, surfaceDir, mLightDir );
        }
        else
        {
            mNormA = new Vector2( 1, 0 );
            mNormB = new Vector2( 1, 0 );
            mNormC = new Vector2( -1, 0 );
            mNormD = new Vector2( -1, 0 );
            mNormE = new Vector2( 0, 1 );
        }

        Vector2.Scale( tmp1, mSurfaceVec, 0.5F );
        Vector2.Add( tmp2, tmp1, mSurfacePoint );
        Vector2.Add( tmp1, mMaxDir, mMinDir );
        Vector2.Scale( tmp1, tmp1, 50 );  //TODO check aliasing
        Vector2.Add( tmp2, tmp1, tmp2 );
        Vector2 core = new Vector2( tmp2 );
        mCoreRgb = ComputeFull( core );        
        System.out.println( "Core = " + core + " " + mCoreRgb );
    }

    int Compute( Vector2 point )
    {
        Vector2.Minus( mTmp1, point, mSurfacePoint );
        if( mTmp1.Dot( mNormE ) >= -1 )
        {
            Vector2.Minus( mTmp2, mTmp1, mSurfaceVec );

            if( mTmp1.Dot( mNormB ) >= 1 && mTmp2.Dot( mNormC ) >= 1 )
            {
                return mCoreRgb;
            }
            else if( mTmp1.Dot( mNormA ) >= -1 && mTmp2.Dot( mNormD ) >= -1 )
            {
                return ComputeFull( point );
            }
            else 
            {
                return 0x000000;
            }
        }
        else
        {
            return 0x000000;
        }
    }

    private int ComputeFull( Vector2 point )
    {
        Vector2.Minus( mTmp1, point, mSurfacePoint );
        Vector2.Scale( mTmp2, mNormB, -0.5F );
        Vector2.Add( mTmp1, mTmp1, mTmp2 );
        float first = (float)Math.min( 1.0, Math.max(
                                           mSurfaceVec.ComputeIntersectScale( mTmp1, mNegB ),
                                           0.0 ) );
        Vector2.Minus( mTmp1, point, mSurfacePoint );
        Vector2.Scale( mTmp2, mNormC, -0.5F );
        Vector2.Add( mTmp1, mTmp1, mTmp2 );
        float last = (float)Math.max( 0.0, Math.min(
                                          mSurfaceVec.ComputeIntersectScale( mTmp1, mNegC ),
                                          1.0 ) );

        int intvls = 10;
        float dist = (last - first) * mSurfaceVec.Norm();
        float weight = dist / intvls;
        float mult = (last - first) / intvls;
        CieXyz xyz = new CieXyz();
        for( int i = 0; i < intvls; ++i )
        {
            Vector2.Scale( mTmp1, mSurfaceVec, first + (i+0.5F)*mult );
            Vector2.Add( mTmp1, mSurfacePoint, mTmp1 );
            Vector2.Minus( mTmp2, point, mTmp1 );
            float distance = mTmp2.Norm();
            float theta = (float)Math.acos(mTmp2.Dot(mNormE)/distance);
            float theta0 = theta - (float)Math.asin(0.5/distance);
            float theta1 = theta + (float)Math.asin(0.5/distance);
            if( (theta0 < mMaxTheta || theta1 > mMinTheta) &&
                mTmp2.Dot( mSurfaceVec ) > 0 )
            {
                float n0 = (float)Math.max(MinN,(Math.sin(theta0)/Math.sin(mThetaD)));
                float n1 = (float)Math.min(MaxN,(Math.sin(theta1)/Math.sin(mThetaD)));
                        
                float w0 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n0 - 4834000,0)));
                float w1 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n1 - 4834000,0)));
                        
                xyz = xyz.AddMult( mObserver.Integrate( w1, w0 ), weight );
            }
        }

        int rgb = xyz.toRgb( 0x0 );
        int r = ((rgb >> 16)%256);
        int g = ((rgb >> 8)%256);
        int b = (rgb %256);
        return (((int)(r*mOrigWeight))<<16) +
            (((int)(g*mOrigWeight))<<8) +
            ((int)(b*mOrigWeight));
    }
/*
                
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
                //rgb = xyz.toRgb( 0x0 );
                
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
				* /
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
*/
}
