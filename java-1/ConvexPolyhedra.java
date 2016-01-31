class ConvexPolyhedra implements IDrawable
{
    static public class Edge
    {
        int a;
        int b;

        Edge( int a_, int b_ )
        {
            a = a_;
            b = b_;
        }
    }

    static public class Face
    {
        int[] edges;
        boolean reverse0;
        boolean reverse1;

        public Face( int a, boolean ra, int b, boolean rb, int c )
        {
            edges = new int[3];
            edges[0] = a;
            edges[1] = b;
            edges[2] = c;
            
            reverse0 = ra;
            reverse1 = rb;
        }

        public Face( int a, boolean ra, int b, boolean rb, int c, int d )
        {
            edges = new int[4];
            edges[0] = a;
            edges[1] = b;
            edges[2] = c;
            edges[3] = d;
            
            reverse0 = ra;
            reverse1 = rb;
        }

        public Face( int a, boolean ra, int b, boolean rb, int c, int d,
                     int e, int f )
        {
            edges = new int[6];
            edges[0] = a;
            edges[1] = b;
            edges[2] = c;
            edges[3] = d;
            edges[4] = e;
            edges[5] = f;
            
            reverse0 = ra;
            reverse1 = rb;
        }

        public Face( int a, boolean ra, int b, boolean rb, int[] r )
        {
            edges = new int[2+r.length];
            edges[0] = a;
            edges[1] = b;
            for( int i = 0; i < r.length; ++i )
            {
                edges[i+2] = r[i];
            }
            
            reverse0 = ra;
            reverse1 = rb;
        }

        public Face( int[] r, boolean ra, boolean rb )
        {
            edges = new int[r.length];
            for( int i = 0; i < r.length; ++i )
            {
                edges[i] = r[i];
            }
            
            reverse0 = ra;
            reverse1 = rb;
        }
    }

    private Vector3[] mPoints;
    private Edge[] mEdges;
    private Face[] mFaces;
    private boolean[] bwork;
    private Vector2 work0;
    private Vector2 work1;

    static public class Parameters
    {
        Vector3[] points;
        Edge[] edges;
        Face[] faces;
    }

    ConvexPolyhedra( Parameters params )
    {
        mPoints = params.points;
        mEdges = params.edges;
        mFaces = params.faces;
        bwork = new boolean[mEdges.length];
        work0 = new Vector2();
        work1 = new Vector2();
    }

    public void Transform( Matrix3 M )
    {
        for( int i = 0; i < mPoints.length; ++i )
        {
            Matrix3.Multiply( mPoints[i], M, mPoints[i] );
        }
    }

    public float ComputeVolume()
    {
        double volume = 0.0;

        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;

        for( int p = 0; p < mPoints.length; ++p )
        {
            Vector3 point = mPoints[p];
            cx += point.x;
            cy += point.y;
            cz += point.z;
        }
        cx /= mPoints.length;
        cy /= mPoints.length;
        cz /= mPoints.length;
        Vector3 o = new Vector3( (float)cx, (float)cy, (float)cz );

        boolean[] reverse = new boolean[8];
        Vector3 e0 = new Vector3();
        Vector3 e1 = new Vector3();
        Vector3 h = new Vector3();
        Vector3 n = new Vector3();

        for( int f = 0; f < mFaces.length; ++f )
        {
            assert mFaces[f].edges.length > 2;

            int sides = mFaces[f].edges.length;
            GetReverse( reverse, f );
            for( int s = 0; s < sides-2; ++s )
            {
                Vector3 a = GetVertex( f, s+0, reverse[s+0] );
                Vector3 b = GetVertex( f, s+1, reverse[s+1] );
                Vector3 c = GetVertex( f, s+2, reverse[s+2] );
                Vector3.Minus( e0, b, a );
                Vector3.Minus( e1, c, b );
                Vector3.Cross( n, e0, e1 );
                Vector3.Minus( h, a, o );
                volume += (1.0/6.0) * h.Dot( n );
            }
        }

        return (float) volume;
    }

    public void Rasterize( float scale, Vector2 transform,
                           RawBitmap bitmap, int rgb )
    {
        for( int i = 0; i < mEdges.length; ++i )
        {
            bwork[i] = false;
        }
        for( int i = 0; i < mFaces.length; ++i )
        {
            if( IsVisible(i) )
            {
                for( int j = 0; j < mFaces[i].edges.length; ++j )
                {
                    bwork[mFaces[i].edges[j]] = true;
                }
            }
        }
        for( int i = 0; i < mEdges.length; ++i )
        {
            if( bwork[i] )
            {
                work0.x = mPoints[mEdges[i].a].x * scale + transform.x;
                work0.y = mPoints[mEdges[i].a].y * scale + transform.y;
                work1.x = mPoints[mEdges[i].b].x * scale + transform.x;
                work1.y = mPoints[mEdges[i].b].y * scale + transform.y;
                Line l = new Line( work0, work1, IDrawable.LineWidth );
                l.Rasterize( bitmap, rgb );
            }
        }
    }

    protected void GetPoint( int index, Vector3 retval )
    {
        Vector3.Copy( retval, mPoints[index] );
    }

    protected void Check()
    {
        System.out.println( "Checking for nulls." );
        for( int i = 0; i < mPoints.length; ++i )
        {
            if( null == mPoints[i] )
            {
                System.out.println( "Point " + i + " is null." );
                System.exit(0);
            }
        }
        for( int i = 0; i < mEdges.length; ++i )
        {
            if( null == mEdges[i] )
            {
                System.out.println( "Edge " + i + " is null." );
                System.exit(0);
            }
        }
        for( int i = 0; i < mFaces.length; ++i )
        {
            if( null == mFaces[i] )
            {
                System.out.println( "Face " + i + " is null." );
                System.exit(0);
            }
        }

        System.out.println( "Checking for degenerate faces." );
        for( int i = 0; i < mFaces.length; ++i )
        {
            if( mFaces[i].edges.length < 3 ) Fail( i );
        }

        System.out.println( "Checking for discontinuous faces." );
        for( int i = 0; i < mFaces.length; ++i )
        {
            boolean[] reverse = GetReverse( i );
            for( int j = 0; j < mFaces[i].edges.length; ++j )
            {
                int jp1 = ((j+1)%mFaces[i].edges.length);
                int end0 = reverse[j] ? mEdges[mFaces[i].edges[j]].a
                    : mEdges[mFaces[i].edges[j]].b;
                int start1 = reverse[jp1] ? mEdges[mFaces[i].edges[jp1]].b
                    : mEdges[mFaces[i].edges[jp1]].a;
                if( end0 != start1 ) Fail( i, j );
            }
        }

        System.out.print( "Checking face colinearity:" );
        float clmx = 0.0F;
        int clf = -1;
        for( int i = 0; i < mFaces.length; ++i )
        {
            Vector3 a = GetVertex( i, 0 );
            Vector3 norm = GetNorm( i );

            for( int j = 0; j < mFaces[i].edges.length - 1; ++j )
            {
                for( int k = 0; k < 2; ++k )
                {
                    int idx = k==0 ? mEdges[mFaces[i].edges[j]].a
                        : mEdges[mFaces[i].edges[j]].b;
                    Vector3 point = mPoints[idx];
                    Vector3 rpoint = new Vector3();
                    Vector3.Minus( rpoint, point, a );
                    Vector3 result = new Vector3();
                    Vector3.IntersectPlane( result, norm, rpoint, norm );
                    Vector3.Minus( result, result, rpoint );
                    if( result.Norm() > clmx )
                    {
                        clmx = result.Norm();
                        clf = i;
                    }
                }
            }
        }
        System.out.println( " " + clmx + " from Face" + clf );

        System.out.print( "Checking convexivity:" );
        float cvmx = 0.0F;
        int cvf = -1;
        for( int i = 0; i < mFaces.length; ++i )
        {
            Vector3 a = GetVertex( i, 0 );
            Vector3 norm = GetNorm( i );

            for( int k = 0; k < mPoints.length; ++k )
            {
                Vector3 rel = new Vector3();
                Vector3.Minus( rel, mPoints[k], a );
                Vector3.Normalize( rel );
                float d = rel.Dot( norm );
                if( d > cvmx )
                {
                    cvmx = d;
                    cvf = i;
                }
            }
        }
        System.out.println( " " + cvmx + " from Face" + cvf );

        System.out.println( "Checking edge inclusion." );
        int[] inc = new int[mEdges.length];
        int[] incr = new int[mEdges.length];
        for( int i = 0; i < mFaces.length; ++i )
        {
            boolean[] reverse = GetReverse( i );
            for( int j = 0; j < mFaces[i].edges.length; ++j )
            {
                int idx = mFaces[i].edges[j];
                if( reverse[j] ) ++incr[idx];
                else ++inc[idx];
            }
        }
        boolean eifail = false;
        for( int i = 0; i < mEdges.length; ++i )
        {
            if( 1 != inc[i] || 1 != incr[i] )
            {
                System.out.println( "Edge" + i + " has face count of " +
                                    (inc[i]+incr[i]) );
                eifail = true;
            }
        }
        if( eifail )
        {
            System.exit(0);
        }

        System.out.println( "Checking Euler Formula." );
        if( mEdges.length != mFaces.length + mPoints.length - 2 )
        {
            System.out.println( "E=" + mEdges.length );
            System.out.println( "F=" + mFaces.length );
            System.out.println( "V=" + mPoints.length );
            System.exit(0);
        }

        for( int i = 0; i < mFaces.length; ++i ){
            System.out.print( "F" + i );
            for( int j = 0; j < mFaces[i].edges.length; ++j ){
                System.out.print( " " + mFaces[i].edges[j] );
            }
            System.out.println();
        }

        for( int i = 0; i < mEdges.length; ++i ){
            System.out.println( "E" + i + " " + mEdges[i].a + " " + mEdges[i].b );
        }

        for( int i = 0; i < mPoints.length; ++i ){
            System.out.println( "P" + i + " " + mPoints[i] );
        }
    }

    private boolean IsVisible( int f )
    {
        Vector3 p0 = mFaces[f].reverse0 
            ? mPoints[mEdges[mFaces[f].edges[0]].b]
            : mPoints[mEdges[mFaces[f].edges[0]].a];
        Vector3 p1 = mFaces[f].reverse0 
            ? mPoints[mEdges[mFaces[f].edges[0]].a]
            : mPoints[mEdges[mFaces[f].edges[0]].b];
        Vector3 p2 = mFaces[f].reverse1 
            ? mPoints[mEdges[mFaces[f].edges[1]].a]
            : mPoints[mEdges[mFaces[f].edges[1]].b];
        float px = p1.x - p0.x;
        float py = p1.y - p0.y;
        float qx = p2.x - p1.x;
        float qy = p2.y - p1.y;
        return px * qy - py * qx > 0;
    }

    private void GetReverse( boolean[] reverse, int face )
    {
        assert reverse.length >= mFaces[face].edges.length;
        reverse[0] = mFaces[face].reverse0;
        reverse[1] = mFaces[face].reverse1;
        for( int i = 2; i < mFaces[face].edges.length; ++i )
        {
            int start = (reverse[i-1]?mEdges[mFaces[face].edges[i-1]].a
                         : mEdges[mFaces[face].edges[i-1]].b);
            reverse[i] = (start != mEdges[mFaces[face].edges[i]].a);
        }
    }

    private boolean[] GetReverse( int face )
    {
        boolean[] retval = new boolean[mFaces[face].edges.length];
        GetReverse( retval, face );
        return retval;
    }

    private Vector3 GetVertex( int face, int point, boolean reverse )
    {
        int idx = reverse ? mEdges[mFaces[face].edges[point]].b
            : mEdges[mFaces[face].edges[point]].a;
        return mPoints[idx];
    }

    private Vector3 GetVertex( int face, int point )
    {
        boolean[] reverse = new boolean[mFaces[face].edges.length];
        GetReverse( reverse, face );
        return GetVertex( face, point, reverse[point] );
    }

    private void Fail( int face )
    {
        String msg = "Face" + face + "(";
        for( int i = 0; i < mFaces[face].edges.length; ++i )
        {
            if( i > 0 )
            {
                msg += ",";
            }
            
            msg += i;
            if( 0 == i && mFaces[face].reverse0 ||
                1 == i && mFaces[face].reverse1 )
            {
                msg += "R";
            }
        }
        msg += ")";
        System.out.println( msg );
        System.exit(-1);
    }

    private void Fail( int face, int edge )
    {
        String msg = "Face" + face + ":Edge" + edge + "(" +
            mEdges[mFaces[face].edges[edge]].a + "," + 
            mEdges[mFaces[face].edges[edge]].b + ")";
        System.out.println( msg );
        System.exit(-1);
    }

    private Vector3 GetNorm( int face )
    {
        Vector3 a = GetVertex( face, 0 );
        Vector3 b = GetVertex( face, 1 );
        Vector3 c = GetVertex( face, 2 );
        Vector3 normal = new Vector3();
        Vector3.Normal( normal, a, b, c );
        return normal;
    }

    static public ConvexPolyhedra MakeTestCube()
    {
        Vector3[] p = new Vector3[8];
        p[0] = new Vector3(0,0,0);
        p[1] = new Vector3(0,0,1);
        p[2] = new Vector3(0,1,0);
        p[3] = new Vector3(0,1,1);
        p[4] = new Vector3(1,0,0);
        p[5] = new Vector3(1,0,1);
        p[6] = new Vector3(1,1,0);
        p[7] = new Vector3(1,1,1);
        Edge e01 = new Edge( 0, 1 );
        Edge e02 = new Edge( 0, 2 );
        Edge e04 = new Edge( 0, 4 );
        Edge e13 = new Edge( 1, 3 );
        Edge e15 = new Edge( 1, 5 );
        Edge e23 = new Edge( 2, 3 );
        Edge e26 = new Edge( 2, 6 );
        Edge e37 = new Edge( 3, 7 );
        Edge e45 = new Edge( 4, 5 );
        Edge e46 = new Edge( 4, 6 );
        Edge e57 = new Edge( 5, 7 );
        Edge e67 = new Edge( 6, 7 );
        Edge[] edges = new Edge[12];
        edges[0]  = e01;
        edges[1]  = e02;
        edges[2]  = e04;
        edges[3]  = e13;
        edges[4]  = e15;
        edges[5]  = e23;
        edges[6]  = e26;
        edges[7]  = e37;
        edges[8]  = e45;
        edges[9]  = e46;
        edges[10] = e57;
        edges[11] = e67;
        Face[] faces = new Face[6];
        faces[0] = new Face( 4, false, 10, false,  7,  3 );
        faces[1] = new Face( 8,  true,  9, false, 11, 10 );
        faces[2] = new Face( 2,  true,  1, false,  6,  9 );
        faces[3] = new Face( 0, false,  3, false,  5,  1 );
        faces[4] = new Face( 5, false,  7, false, 11,  6 );
        faces[5] = new Face( 2, false,  8, false,  4,  0 );

        Parameters params = new Parameters();
        params.points = p;
        params.edges = edges;
        params.faces = faces;
        return new ConvexPolyhedra( params );
    }
}
