public class Prism implements PixelRasterizer
{
    private enum Region{ I, II, III, IV };
    static private final float Ninety = (float)(Math.PI / 2.0);
    static private final double ExSinTheta = Math.cos( 0.0001 );

    public final Vector2 mSurfacePoint;
    public final Vector2 mSurfaceVec;
    public final Vector2 mMinDir;
    public final Vector2 mMaxDir;

    private final Material mMaterial;
    private final CieXyzObserver mObserver;
    private final float mWeight;
    private final float mThetaD;
    private final float mMinTheta;
    private final float mMaxTheta;

    private final Vector2 mNegB;
    private final Vector2 mNegC;

    private final Vector2 mNormA;
    private final Vector2 mNormB;
    private final Vector2 mNormC;
    private final Vector2 mNormD;
    private final Vector2 mNormE;

    private final boolean mHasCoreSample;
    private final int mCoreRgb;

    private Vector2 mTmp1;
    private Vector2 mTmp2;
    private Vector2 mTmp3;

    public Prism( Material material,
                  CieXyzObserver observer,
                  Vector2 surface1,
                  Vector2 surface2,
                  Vector2 lightDir,
                  float weight )
    {
        mMaterial = material;
        mObserver = observer;
        mWeight = weight;

        Vector2 tmp1 = new Vector2();
        Vector2 tmp2 = new Vector2();

        Vector2.Minus( tmp1, surface2, surface1 );
        if( lightDir.Dot( tmp1 ) >= 0.0F )
        {
            mSurfacePoint = new Vector2( surface1 );
            Vector2.Minus( tmp1, surface2, surface1 );
            mSurfaceVec = new Vector2(tmp1);
        }
        else
        {
            mSurfacePoint = new Vector2( surface2 );
            Vector2.Minus( tmp1, surface1, surface2 );
            mSurfaceVec = new Vector2(tmp1);
        }

        Vector2 surfaceDir;
        if( mSurfaceVec.Norm2() > 1e-8F )
        {
            surfaceDir = new Vector2( mSurfaceVec );
            Vector2.Normalize( surfaceDir );
        }
        else
        {
            surfaceDir = new Vector2( lightDir );
            Vector2.Normalize( surfaceDir );
        }
        mThetaD = Ninety - 
            (float) Math.acos( surfaceDir.Dot( lightDir ) /
                               lightDir.Norm());
        Vector2 normal = new Vector2( surfaceDir.y, -surfaceDir.x );
        Vector2.Normalize( normal );
        if( lightDir.Dot( normal ) < 0 )
        {
            Vector2.Scale( normal, normal, -1 );
        }

        if( surfaceDir.Dot( lightDir ) > 1E-8 )
        {
            Vector2.NormFromRel( tmp2, surfaceDir, lightDir );
            mNormE = new Vector2( tmp2 );

            double minSinTheta = mMaterial.MinN*Math.sin(mThetaD);
            if( -ExSinTheta < minSinTheta && minSinTheta < ExSinTheta )
            {
                mMinTheta = (float)Math.asin(mMaterial.MinN*Math.sin(mThetaD));
                Vector2.Scale( tmp1, surfaceDir, (float)Math.sin( mMinTheta ) );
                Vector2.Scale( tmp2, normal, (float)Math.cos( mMinTheta ) );
                Vector2.Add( tmp1, tmp1, tmp2 );
                Vector2.Normalize( tmp1 );
                mMinDir = new Vector2( tmp1 );

                Vector2.Scale( tmp1, surfaceDir, -1 );
                Vector2.NormFromRel( tmp2, mMinDir, surfaceDir );
                mNormA = new Vector2( tmp2 );
                Vector2.NormFromRel( tmp2, mMinDir, tmp1 );
                mNormC = new Vector2( tmp2 );
            }
            else
            {
                mMinTheta = Ninety;
                mMinDir = new Vector2( surfaceDir );
                Vector2.Scale( tmp1, lightDir, -1 );
                Vector2.NormFromRel( tmp2, mMinDir, tmp1 );
                mNormA = new Vector2( tmp2 );
                Vector2.NormFromRel( tmp2, mMinDir, lightDir );
                mNormC = new Vector2( tmp2 );
            }
            Vector2.Scale( tmp2, mMinDir, -1/mMinDir.Norm() );
            mNegC = new Vector2( tmp2 );

            double maxSinTheta = mMaterial.MaxN*Math.sin(mThetaD);
            if( -ExSinTheta < maxSinTheta && maxSinTheta < ExSinTheta )
            {
                mMaxTheta = (float)Math.asin(mMaterial.MaxN*Math.sin(mThetaD));
                Vector2.Scale( tmp1, surfaceDir, (float)Math.sin( mMaxTheta ) );
                Vector2.Scale( tmp2, normal, (float)Math.cos( mMaxTheta ) );
                Vector2.Add( tmp1, tmp1, tmp2 );
                Vector2.Normalize( tmp1 );
                mMaxDir = new Vector2( tmp1 );
 
                Vector2.Scale( tmp1, surfaceDir, -1 );
                Vector2.NormFromRel( tmp2, mMaxDir, surfaceDir );
                mNormB = new Vector2( tmp2 );
                Vector2.NormFromRel( tmp2, mMaxDir, tmp1 );
                mNormD = new Vector2( tmp2 );
            }
            else
            {
                mMaxTheta = Ninety;
                mMaxDir = new Vector2( surfaceDir );

                Vector2.Scale( tmp1, lightDir, -1 );
                Vector2.NormFromRel( tmp2, mMaxDir, tmp1 );
                mNormB = new Vector2( tmp2 );
                Vector2.NormFromRel( tmp2, mMaxDir, lightDir );
                mNormD = new Vector2( tmp2 );
            }
            Vector2.Scale( tmp2, mMaxDir, -1/mMaxDir.Norm() );
            mNegB = new Vector2( tmp2 );
        }
        else
        {
            mMinTheta = Ninety;
            mMaxTheta = Ninety;
            mMinDir = new Vector2( 1, 0 );
            mMaxDir = new Vector2( 1, 0 );
            mNormA = new Vector2( 1, 0 );
            mNormB = new Vector2( 1, 0 );
            mNormC = new Vector2( -1, 0 );
            mNormD = new Vector2( -1, 0 );
            mNormE = new Vector2( 0, 1 );
            mNegB = new Vector2( -1, 0 );
            mNegC = new Vector2( -1, 0 );
        }

        mTmp1 = new Vector2();
        mTmp2 = new Vector2();
        mTmp3 = new Vector2();

        float region3Tip = mMaxDir.ComputeIntersectScale( mSurfaceVec,
                                                          mMinDir );
        float coreSampleLen;
        if( Math.abs( region3Tip ) < 50 )
        {
            coreSampleLen = (float) Math.abs( region3Tip ) / 2.0F;
        }
        else
        {
            coreSampleLen = 50;
        }

        Vector2.Scale( tmp1, mSurfaceVec, 0.5F );
        Vector2.Add( tmp1, tmp1, mSurfacePoint );
        Vector2.Add( tmp2, mMaxDir, mMinDir );
        Vector2.Scale( tmp2, tmp2, coreSampleLen );        
        Vector2.Add( tmp1, tmp1, tmp2 );
        if( Region.III == GetRegion( tmp1 ) )
        {
            mHasCoreSample = true;
            mCoreRgb = ComputeFull( tmp1 );
        }
        else
        {
            mHasCoreSample = false;
            mCoreRgb = 0;
        }
    }

    public int PixelRasterize( int x, int y )
    {
        Vector2 point = new Vector2( x, y );
        return Compute( point, GetRegion( point ) );
    }

    public ConvexPolygon Intersect( ConvexPolygon bounds )
    {
        if( mMinTheta >= Ninety - 1e-8 )
        {
            assert mMaxTheta >= Ninety - 1e-8;
            return ConvexPolygon.Empty;
        }

        Vector2.Scale( mTmp1, mSurfaceVec, 0.5F );
        Vector2.Add( mTmp1, mSurfacePoint, mTmp1 );
        
        bounds.GetCenter( mTmp2 );
        Vector2.Minus( mTmp2, mTmp2, mTmp1 );
        float radius = bounds.GetRadius() + mTmp2.Norm();

        Vector2.Scale( mTmp1, mSurfaceVec, 0.5F );
        Vector2.Minus( mTmp1, mTmp1, mMinDir );

        Vector2.Add( mTmp2, mSurfaceVec, mMaxDir );
        Vector2.Minus( mTmp2, mTmp2, mMinDir );
        Vector2.Normalize( mTmp2 );

        Vector2.Scale( mTmp2, mTmp2, mTmp1.Dot( mTmp2 ) );
        Vector2.Minus( mTmp1, mTmp1, mTmp2 );
        
        System.out.println( "H = " + mTmp1.Norm() );
        float length = radius / mTmp1.Norm();
        length = (float)Math.min( length, 1000 );

        Vector2.Add( mTmp1, mSurfacePoint, mSurfaceVec );
        Vector2.Scale( mTmp2, mMaxDir, length );
        Vector2.Add( mTmp2, mTmp1, mTmp2 );
        Vector2.Scale( mTmp3, mMinDir, length );
        Vector2.Add( mTmp3, mSurfacePoint, mTmp3 );

        Vector2[] points;
        if( mMaxTheta <= Ninety - 1e-8 )
        {
            points = new Vector2[4];
            if( mSurfaceVec.Cross( mNormE ) > 0 )
            {
                points[0] = mSurfacePoint;
                points[1] = mTmp1;
                points[2] = mTmp2;
                points[3] = mTmp3;
            }
            else
            {
                points[3] = mSurfacePoint;
                points[2] = mTmp1;
                points[1] = mTmp2;
                points[0] = mTmp3;            
            }
        }
        else
        {
            points = new Vector2[3];
            if( mSurfaceVec.Cross( mNormE ) > 0 )
            {
                points[0] = mSurfacePoint;
                points[1] = mTmp2;
                points[2] = mTmp3;
            }
            else
            {
                points[2] = mSurfacePoint;
                points[2] = mTmp2;
                points[0] = mTmp3;            
            }
        }

        return ConvexPolygon.Intersect( new ConvexPolygon( points ), bounds );
    }

    private Region GetRegion( Vector2 point )
    {
        Vector2.Minus( mTmp1, point, mSurfacePoint );
        if( mTmp1.Dot( mNormE ) >= -1 )
        {
            Vector2.Minus( mTmp2, mTmp1, mSurfaceVec );

            if( mTmp1.Dot( mNormA ) < -1 || mTmp2.Dot( mNormD ) < -1 )
            {
                return Region.II;
            }
            else if( mTmp1.Dot( mNormB ) >= 1 && mTmp2.Dot( mNormC ) >= 1 )
            {
                return Region.III;
            }
            {
                return Region.IV;
            }
        }
        else
        {
            return Region.I;
        }
    }

    private int Compute( Vector2 point, Region region )
    {
        switch( region )
        {
        case I:
        case II:
            return 0x000000;
        case III:
            if( mHasCoreSample )
            {
                return mCoreRgb;
            }
            else
            {
                return ComputeFull( point );
            }
        case IV:
            return ComputeFull( point );
        }
        return 0x000000;
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
                float n0 = (float)Math.max(mMaterial.MinN,(Math.sin(theta0)/Math.sin(mThetaD)));
                float n1 = (float)Math.min(mMaterial.MaxN,(Math.sin(theta1)/Math.sin(mThetaD)));
                        
                float w0 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n0 - 4834000,0)));
                float w1 = (float)(731.3 - Math.sqrt(Math.max(2008000 * n1 - 4834000,0)));
                        
                xyz = xyz.AddMult( mObserver.Integrate( w1, w0 ), weight );
            }
        }

        int rgb = xyz.toRgb( 0x0 );
        int r = ((rgb >> 16)%256);
        int g = ((rgb >> 8)%256);
        int b = (rgb %256);
        return (((int)(r*mWeight))<<16) +
            (((int)(g*mWeight))<<8) +
            ((int)(b*mWeight));
    }
}
