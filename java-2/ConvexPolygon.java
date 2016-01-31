import gameoxide.GOImageCanvas;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

class ConvexPolygon
{
    public ConvexPolygon( Vector2[] points, Material m )
    {
        mMaterial = new Material( m );
        mPoints = new Vector2[points.length];
        mVecs = new Vector2[points.length];
        mNorms = new Vector2[points.length];
        for( int i = 0; i < points.length; ++i )
        {
            int next = (i+1)%points.length;
            mPoints[i] = new Vector2( points[i] );
            mVecs[i] = new Vector2();
            Vector2.Minus( mVecs[i], points[next], points[i] );
        }
        for( int i = 0; i < points.length; ++i )
        {
            int next = (i+1)%points.length;
            mNorms[i] = new Vector2();
            Vector2.NormFromRel( mNorms[i], mVecs[i], mVecs[next] );
            System.out.println( mVecs[i] + " and " + mVecs[next] + " norm " + mNorms[i] );
        }

        mTmp1 = new Vector2();
        mTmp2 = new Vector2();
        mTmp3 = new Vector2();
    }

    List<LightBeam> Process( LightBeam beam )
    {
        System.out.println("Processing " + beam);
        List<LightBeam> retval = new LinkedList<LightBeam>();
        for( int i = 0; i < mPoints.length; ++i )
        {
            System.out.print("  Checking against " + mVecs[i] );
            if( beam.dir.Dot( mNorms[i] ) < -1E-4F )
            {
                Vector2.Minus( mTmp1, beam.src.start, mPoints[i] );
                float s = mVecs[i].ComputeIntersectScale( mTmp1, beam.dir );
                Vector2.Minus( mTmp1, beam.src.end, mPoints[i] );
                float t = mVecs[i].ComputeIntersectScale( mTmp1, beam.dir );
                
                //Vector2.Minus( mTmp1, beam.src.end, beam.src.start );
                //if( mTmp1.Dot( mVecs[i] ) < 0 )
                if( s > t )
                {
                    float k = s;
                    s = t;
                    t = k;
                }

                if( s < 1 - 1E-4F && t > 1E-4F )
                {
                    s = (float)Math.max(s,0);
                    t = (float)Math.min(t,1);

                    LightBeam lb = new LightBeam();
                    lb.dir = new Vector2();
                    lb.weight = beam.weight;
                    lb.src = new LightBeam.CrossSection();
                    lb.src.start = new Vector2();
                    lb.src.end = new Vector2();
                    lb.src.weight = beam.weight;
                    Vector2.Scale( mTmp1, mVecs[i], s );
                    Vector2.Add( lb.src.start, mPoints[i], mTmp1 );
                    Vector2.Scale( mTmp1, mVecs[i], t );
                    Vector2.Add( lb.src.end, mPoints[i], mTmp1 );

                    Vector2.Minus( mTmp1, beam.src.end, beam.src.start );
                    Vector2.Minus( mTmp2, lb.src.start, beam.src.start );
                    Vector2.Scale( mTmp3, beam.dir, -1 );
                    float lbs = mTmp1.ComputeIntersectScale( mTmp2, mTmp3 ); 
                    Vector2.Minus( mTmp2, lb.src.end, beam.src.start );
                    float lbt = mTmp1.ComputeIntersectScale( mTmp2, mTmp3 ); 
                    lb.path = beam.scalePath( lbs, lbt, mTmp1 );

                    float cosTheta = -mNorms[i].Dot( beam.dir );
                    float reflectedWeight;
                    float transmittedWeight;
                    if( cosTheta <= mMaterial.CosTheta )
                    {
                        reflectedWeight = beam.weight;
                        transmittedWeight = 0;
                    }
                    else
                    {
                        reflectedWeight = mMaterial.EstimateReflectionWeight(
                            (float)Math.acos( cosTheta ) ) * beam.weight;
                        transmittedWeight = beam.weight - reflectedWeight;
                    }

                    if( transmittedWeight > 0 )
                    {
                        LightBeam tlb = new LightBeam();
                        tlb.src = new LightBeam.CrossSection();
                        tlb.src.start = new Vector2( lb.src.start );
                        tlb.src.end = new Vector2( lb.src.end );
                        tlb.path = lb.path;
                        lb.path = new LinkedList<LightBeam.CrossSection>();
                        tlb.weight = transmittedWeight;
                        tlb.inside = false;
                        tlb.dir = beam.dir;
                        retval.add( tlb );
                    }

                    if( reflectedWeight > 0 )
                    {
                        Vector2.Scale( mTmp1, mNorms[i], 
                                       beam.dir.Dot(mNorms[i]) );
                        Vector2.Minus( mTmp2, beam.dir, mTmp1 );
                        Vector2.Minus( lb.dir, mTmp2, mTmp1 );
                        Vector2.Normalize( lb.dir );
                        lb.weight = reflectedWeight;
                        lb.inside = true;

                        retval.add( lb );
                    }
                    System.out.println( " intersected." );
                }
                else
                {
                    System.out.println( " outside bounds (" + s + "," + t + ")" );
                }
            }
            else
            {
                System.out.println( " wrong side." );
            }
        }
        return retval;
    }

    void Draw( GOImageCanvas canvas )
    {
        canvas.setColor( 0x00FFFF );
        for( int i = 0; i < mPoints.length; ++i )
        {
            int next = (i+1)%mPoints.length;
            canvas.moveTo( (int)mPoints[i].x, (int)mPoints[i].y );
            canvas.drawTo( (int)mPoints[next].x, (int)mPoints[next].y );
        }
    }

    Material mMaterial;
    Vector2[] mPoints;
    Vector2[] mVecs;
    Vector2[] mNorms;
    Vector2 mTmp1;
    Vector2 mTmp2;
    Vector2 mTmp3;
}