import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

class ConvexPolygon
{
    static public final ConvexPolygon Empty = new ConvexPolygon();

    private Vector2[] mPoints;
    private Vector2[] mVecs;
    private Vector2 mCenter;
    private Float mRadius2;

    public ConvexPolygon( Vector2[] points )
    {
        mPoints = new Vector2[points.length];
        mVecs = new Vector2[points.length];
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
            if( mVecs[i].Cross( mVecs[next] ) <= 0 )
            {
                System.out.println( i + " " + this );
            }
            assert mVecs[i].Cross( mVecs[next] ) > 0; 
        }
    }

    public ConvexPolygon( ConvexPolygon p )
    {
        mPoints = new Vector2[p.mPoints.length];
        mVecs = new Vector2[p.mPoints.length];
        for( int i = 0; i < p.mPoints.length; ++i )
        {
            mPoints[i] = new Vector2( p.mPoints[i] );
            mVecs[i] = new Vector2( p.mVecs[i] );
        }
        if( p.mCenter != null )
        {
            mCenter = new Vector2( p.mCenter );
            mRadius2 = new Float( p.mRadius2 );
        }
    }

    private ConvexPolygon()
    {
        mPoints = new Vector2[0];
        mVecs = new Vector2[0];
    }

    public int GetSides()
    {
        return mPoints.length;
    }

    public void GetPoint( int i, Vector2 r )
    {
        r.x = mPoints[i].x;
        r.y = mPoints[i].y;
    }

    public void GetCenter( Vector2 c )
    {
        if( mCenter == null ) Initialize();
        c.x = mCenter.x;
        c.y = mCenter.y;
    }

    public float GetRadius()
    {
        if( mCenter == null ) Initialize();
        return (float)Math.sqrt( mRadius2 );
    }

    public void Rasterize( RawBitmap bitmap, int rgb )
    {
        for( int i = 0; i < mPoints.length-2; ++i )
        {
            Triangle t = new Triangle( mPoints[0],
                                       mPoints[i+1],
                                       mPoints[i+2] );
            t.Rasterize( bitmap, rgb );
        }
    }

    public void Rasterize( RawBitmap bitmap, IRasterizer pixels )
    {
        int lowestIndex = 0;
        float lowestY = mPoints[0].y;
        int highestIndex = 0;
        float highestY = mPoints[0].y;
        for( int i = 0; i < mPoints.length; ++i )
        {
            if( mPoints[i].y < lowestY )
            {
                lowestIndex = i;
                lowestY = mPoints[i].y;
            }
            if( mPoints[i].y > highestY )
            {
                highestIndex = i;
                highestY = mPoints[i].y;
            }
        }

        System.out.print( "POLYGON" );
        for( int i = 0, _i = lowestIndex; i < mPoints.length; ++i, _i=(_i+1)%mPoints.length )
        {
            System.out.print( " " + mPoints[_i] );
        }
        System.out.println();

        int yStart = (int)(lowestY+0.5F);
        int yEnd = (int)(highestY+0.5F);
        int left = lowestIndex;
        int right = lowestIndex;
        for( int y = yStart; y <= yEnd; ++y )
        {
            float currTop = y + 0.5F;
            float currBottom = y - 0.5F;

            int nextLeft = (left + mPoints.length-1)%mPoints.length;
            float mLeft = (mPoints[nextLeft].x - mPoints[left].x) /
                (mPoints[nextLeft].y - mPoints[left].y);
            float minX = (currBottom - mPoints[left].y) * mLeft+
                mPoints[left].x;
            boolean leftAdvanced = false;

            do
            {
                if( mPoints[nextLeft].y < currTop )
                {
                    left = nextLeft;
                    if( mPoints[left].x < minX ) minX = mPoints[left].x;
                    nextLeft = (left + mPoints.length-1)%mPoints.length;
                    leftAdvanced = true;
                }
                else
                {
                    if( leftAdvanced )
                    {
                        mLeft = (mPoints[nextLeft].x - mPoints[left].x) /
                            (mPoints[nextLeft].y - mPoints[left].y);
                    }
                    float leftTop = (currTop - mPoints[left].y) * mLeft +
                        mPoints[left].x;
                    if( leftTop < minX ) minX = leftTop;
                    break;
                }
            } while( left != highestIndex );

            int nextRight = (right + 1)%mPoints.length;
            float mRight = (mPoints[nextRight].x - mPoints[right].x) /
                (mPoints[nextRight].y - mPoints[right].y);
            float maxX = mPoints[right].x +
                (currBottom - mPoints[right].y) * mRight;
            boolean rightAdvanced = false;

            do
            {
                if( mPoints[nextRight].y < currTop )
                {
                    right = nextRight;
                    if( mPoints[right].x > maxX ) maxX = mPoints[right].x;
                    nextRight = (right + 1)%mPoints.length;
                    rightAdvanced = true;
                }
                else
                {
                    if( rightAdvanced )
                    {
                        mRight = (mPoints[nextRight].x - mPoints[right].x) /
                            (mPoints[nextRight].y - mPoints[right].y);
                    }
                    float rightTop = mPoints[right].x +
                        (currTop - mPoints[right].y) * mRight;
                    if( rightTop > maxX ) maxX = rightTop;
                    break;
                }
            } while( right != highestIndex );

            int xStart = (int)(minX+0.5F);
            int xEnd = (int)(maxX+0.5F);
            for( int x = xStart; x <= xEnd; ++x )
            {
                bitmap.Add( x, y, pixels.GetPixel( x, y ) );
            }
        }
    }

    public String toString()
    {
        String retval = "POLY " + mPoints.length + ":"; 
        for( int i = 0; i < mPoints.length; ++i )
        {
            retval += " " + mPoints[i];
        }
        return retval;
    }

    private void Initialize()
    {
        float xb = 0.0F;
        float yb = 0.0F;
        for( int i = 0; i < mPoints.length; ++i )
        {
            xb = xb + (mPoints[i].x - xb) / (i+1);
            yb = yb + (mPoints[i].y - yb) / (i+1);
        }
        mCenter = new Vector2( xb, yb );

        Vector2 tmp = new Vector2();
        float r2 = 0.0F;
        for( int i = 0; i < mPoints.length; ++i )
        {
            Vector2.Minus( tmp, mPoints[i], mCenter );
            float c2 = tmp.Norm2();
            if( c2 > r2 )
            {
                r2 = c2;
            }
        }
        mRadius2 = new Float( r2 ); 
    }

    static public ConvexPolygon Intersect( ConvexPolygon p, ConvexPolygon q )
    {
        if( q.mPoints.length < p.mPoints.length )
        {
            ConvexPolygon tmp = p;
            p = q;
            q = tmp;
        }

        Vector2 tmp1 = new Vector2();

        int i = 0;
        int j = 0;
        int cnt = 0;

        System.out.println( "START" );
        while( cnt < p.mPoints.length + q.mPoints.length )
        {
            boolean piq = Inside( p, i, q, j, tmp1 );
            boolean qip = Inside( q, j, p, i, tmp1 );

            if( !piq && !qip )
            {
                i = (i+1)%p.mPoints.length;
                j = (j+1)%q.mPoints.length;
            }
            else 
            {
                if( !piq )
                {
                    ConvexPolygon swap = p; p = q; q = swap;
                    int swapi = i; i = j; j = swapi;
                }

                float cp = Scale( p.mPoints[i], p.mVecs[i],
                                  q.mPoints[j], q.mVecs[j], tmp1 );
                float cq = Scale( q.mPoints[j], q.mVecs[j],
                                  p.mPoints[i], p.mVecs[i], tmp1 );

                if( cp >= 0 && cp < 1 && cq >= 0 && cq < 1 )
                {
                    return Iterate(p, i, q, j, cq );
                }
                else
                {
                    if( cp < 1 && cq >= 0 )
                    {
                        System.out.println( " P INSIDE Q, Q OVERLAP" );
                        i = (i+1)%p.mPoints.length;
                    }
                    else
                    {
                        System.out.println( " P INSIDE Q" );
                        j = (j+1)%q.mPoints.length;
                    }
                }
            }

            ++cnt;
        }

        if( p.mCenter == null ) p.Initialize();
        if( q.mCenter == null ) q.Initialize();

        if( q.mRadius2 < p.mRadius2 )
        {
            ConvexPolygon swap = p; p = q; q = swap;
        }

        for( int k = 0; k < q.mPoints.length; ++k )
        {
            Vector2.Minus( tmp1, p.mCenter, q.mPoints[k] );
            if( q.mVecs[k].Cross( tmp1 ) <= 0 )
            {
                return Empty;
            }
        }

        return new ConvexPolygon( p );
    }

    static private ConvexPolygon Iterate( ConvexPolygon p, int i0,
                                          ConvexPolygon q, int j0,
                                          float cqInit )
    {
        Vector2 tmp1 = new Vector2();

        List<Vector2> r = new LinkedList<Vector2>();
        Vector2 init = new Vector2();
        Vector2.Scale( tmp1, p.mVecs[i0], cqInit );
        Vector2.Add( init, p.mPoints[i0], tmp1 );
        r.add( init );

        int icnt = 1;
        int jcnt = 0;
        int i = (i0+1)%p.mPoints.length;
        int j = j0;
        ConvexPolygon swap; int swapi;
        swap = p; p = q; q = swap; 
        swapi = i; i = j; j = swapi;
        swapi = icnt; icnt = jcnt; jcnt = swapi;

        while( icnt < p.mPoints.length )
        {
            if( p.mVecs[i].Cross( q.mVecs[j] ) <= 0 )
            {
                System.out.println( "CROSS j = " + j );
                j = (j+1)%q.mPoints.length;
                ++jcnt;
            }
            else
            {
                float cp = Scale( p.mPoints[i], p.mVecs[i],
                                  q.mPoints[j], q.mVecs[j], tmp1 );
                float cq = Scale( q.mPoints[j], q.mVecs[j],
                                  p.mPoints[i], p.mVecs[i], tmp1 );

                if( cq < 1 && cp < 1)
                {
                    System.out.println( "ISECT i = " + i );
                    assert cq > -1.0E-2F;
                    Vector2 next = new Vector2();
                    Vector2.Scale( next, p.mVecs[i], cq );
                    Vector2.Add( next, p.mPoints[i], next );
                    r.add( next );

                    swap = p; p = q; q = swap; 
                    swapi = i; i = j; j = swapi;
                    swapi = icnt; icnt = jcnt; jcnt = swapi;
                    
                    j = (j+1)%q.mPoints.length;
                    ++jcnt;
                }
                else 
                {
                    if( cp < 1 )
                    {                        
                        System.out.println( "OVER  i = " + i );
                        i = (i+1)%p.mPoints.length;
                        r.add( p.mPoints[i] );
                        ++icnt;
                    }
                    else
                    {
                        System.out.println( "CURL  j = " + j );
                        j = (j+1)%q.mPoints.length;
                        ++jcnt;
                    }
                }
            }                    
        }

        Vector2[] ra = new Vector2[r.size()];
        ListIterator<Vector2> it = r.listIterator();
        int k = 0;
        System.out.println( "Intersection" );
        for( int a = 0; a < p.mPoints.length; ++a )
            System.out.println( p.mPoints[a].toString() );
        System.out.println( " WITH " );
        for( int a = 0; a < q.mPoints.length; ++a )
            System.out.println( q.mPoints[a].toString() );
        System.out.println( " IS " );
        while( it.hasNext() )
        {
            ra[k] = it.next();
            System.out.println( ra[k].toString() );
            ++k;
        }
        return new ConvexPolygon( ra );
    }

    static private boolean Inside( ConvexPolygon p, int i,
                                   ConvexPolygon q, int j, 
                                   Vector2 w )
    {
        Vector2.Minus( w, p.mPoints[i], q.mPoints[j] );
        return q.mVecs[j].Cross( w ) > 0;
    }

    static private float Scale( Vector2 p1, Vector2 v1, Vector2 p2,
                                Vector2 v2, Vector2 w )
    {
        Vector2.Minus( w, p1, p2 );
        return v2.ComputeIntersectScale( w, v1 ); 
    }
}