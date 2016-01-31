class RoundCut extends ConvexPolyhedra implements ICrossSection
{
    static public class Measurements
    {
        Float minWidth;
        Float maxWidth;
        Float height;
        Float tablePct;
        Float heightPct;
        Float crownHeightPct;
        Float pavilionHeightPct;
        Float crownAngleDeg;
        Float pavilionAngleDeg;
        Float girdlePct;
        Float starPct;
        Float pavilionPct;
        Integer girdleEdges;
        boolean hasCulet;

        public void Complete()
        {
            heightPct = new Float(
                height.floatValue() * 2.0F /
                (minWidth.floatValue() + maxWidth.floatValue()) );
            crownHeightPct = new Float(
                0.5F*(1.0F - tablePct.floatValue()) *
                (float)Math.tan( crownAngleDeg.floatValue() * Math.PI/180.0 ));
            pavilionHeightPct = new Float(
                0.5F *
                (float)Math.tan( pavilionAngleDeg.floatValue()*Math.PI/180.0));
            girdlePct = new Float(
                heightPct.floatValue() - 
                crownHeightPct.floatValue() -
                pavilionHeightPct.floatValue() );
        }

        boolean IsComplete()
        {
            return minWidth != null &&
                maxWidth != null &&
                height != null &&
                tablePct != null &&
                heightPct != null &&
                crownHeightPct != null &&
                pavilionHeightPct != null &&
                crownAngleDeg != null &&
                pavilionAngleDeg != null &&
                girdlePct != null &&
                starPct != null &&
                pavilionPct != null &&
                girdleEdges != null;
        }
    }

    public RoundCut( Measurements m )
    {
        super( ComputeParameters(m) );
        mHasCulet = m.hasCulet;
        mIPoints = 4 + 3 + 4*m.girdleEdges + (mHasCulet?1:0);
        mNPoints = 8*mIPoints + (mHasCulet?0:1);
        Check();
    }

    public int GetCrossSectionCnt()
    {
        return 1;//2;
    }

    public Matrix3 GetCrossSectionRot( int index )
    {
        return sRot[index];
    }

    public Vector3[] GetCrossSection( int index )
    {
        Vector3[] points = new Vector3[7 + (mHasCulet?1:0)];
        points[0] = CopyPoint( mHasCulet?mIPoints-1:8*mIPoints );
        points[1] = CopyPoint( PT );
        points[2] = CopyPoint( BB );
        points[3] = CopyPoint( BT );
        points[4] = CopyPoint( BT+mIPoints*4 );
        points[5] = CopyPoint( BB+mIPoints*4 );
        points[6] = CopyPoint( PT+mIPoints*4 );
        if( mHasCulet )
        {
            points[0] = CopyPoint( mIPoints-1+mIPoints*4 );
        }
        return points;
    }

    private boolean mHasCulet;
    private int mNPoints;
    private int mIPoints;

    private Vector3 CopyPoint( int index )
    {
        Vector3 rt = new Vector3();
        GetPoint( index, rt );
        return rt;
    }

    static private Matrix3[] sRot;

    static
    {
        sRot = new Matrix3[2];
        sRot[0] = Matrix3.MakeIdentity();
        sRot[1] = Matrix3.MakeYRotation( (float)(2.0 * Math.PI / 16.0 ) );
    }

    static final int BT = 0;
    static final int BR = 1;
    static final int BB = 2;
    static final int U0R = 3;
    static final int PT = 4;
    static final int PR = 5;
    static final int L0R = 6;

    static private Parameters ComputeParameters( Measurements m )
    {
        /*
         * Naming Convention:
         * 
         * {facet}{position}
         * 
         * Facets:
         * s  : Star facets - Triangle connected to top table.
         * b  : Bezel facets - Kite shapes, upper formed by star facets.
         * u# : Upper girdle - Triangles connected to upper girdle.
         * l# : Lower girdle - Triangles connected to lower girdle.
         * p  : Pavilion - Kite shapes, connected to culet point.
         *
         * Positions:
         * t  : Top point - For use with bezel, upper girdle and pavilion.
         * l  : Left point - For use with all.
         * i# : Internal points - For use with girdle connections. 
         * r  : Right point - For use with all.
         * b  : Bottom point - For use with bezel, lower facets.
         */              

        assert m.IsComplete();
        float tableh = m.tablePct.floatValue()/2.0F;
        float crown = m.crownHeightPct.floatValue();
        float star = m.starPct.floatValue();
        float pavilion = m.pavilionPct.floatValue();
        float girdle = m.girdlePct.floatValue();
        int girdleEdges = m.girdleEdges.intValue();
        float depth = girdle + m.pavilionHeightPct.floatValue();
        float thetad = (float)(2.0 * Math.PI / 8.0 );
        boolean hasCulet = false;

        final int ipoints = 4 + 3 + 4*girdleEdges + (hasCulet?1:0);
        final int npoints = 8*ipoints + (hasCulet?0:1);
        Vector3[] points = new Vector3[npoints];
        final int iedges = 6 + 4 + 4*(girdleEdges+1) + 2*girdleEdges + (hasCulet?1:0);
        final int nedges = 8*iedges;
        Edge[] edges = new Edge[nedges];
        final int ifaces = 4 + 3 + 2*girdleEdges;
        final int nfaces = 8*ifaces + 1 + (hasCulet?1:0);
        Face[] faces = new Face[nfaces];
        // E = 8*(6+4+4(g+1)+2g+c) = 112+48g+8c
        // F = 8*(4+3+2g)+1+c = 57+16g+c
        // V = 8*(4+3+4g+c)+(1-c) = 57+32g+7c
        // F+V-2 = 112+48g+8c

        Vector3 lp[] = new Vector3[ipoints];

        final int U0Iref = 7;
        final int U1Iref = U0Iref + girdleEdges;
        final int L0Iref = U1Iref + girdleEdges;
        final int L1Iref = L0Iref + girdleEdges;
        final int CL = (hasCulet ? ipoints-1 : -1);

        //Compute primary vectors
        Vector2 a0 = new Vector2( 1.0F, 0.0F );
        Vector2 a1 = new Vector2( (float)Math.cos(thetad/2.0),
                                  (float)Math.sin(thetad/2.0) );
        Vector2 a2 = new Vector2( (float)Math.cos(thetad),
                                  (float)Math.sin(thetad) );

        //Compute major markers.
        lp[BT] = new Vector3( tableh*a0.x, crown, tableh*a0.y );
        lp[BB] = new Vector3( 0.5F*a0.x, 0.0F, 0.5F*a0.y );
        lp[U0R] = new Vector3( 0.5F*a1.x, 0.0F, -0.5F*a1.y );
        Vector3 u1r = new Vector3( 0.5F*a2.x, 0.0F, -0.5F*a2.y );
        lp[PT] = new Vector3( 0.5F*a0.x, -girdle, 0.5F*a0.y );
        lp[L0R] = new Vector3( 0.5F*a1.x, -girdle, -0.5F*a1.y );
        Vector3 l1r = new Vector3( 0.5F*a2.x, -girdle, -0.5F*a2.y );
        Vector3 pb = new Vector3( 0.0F, -depth, 0.0F );

        //Compute intersect vector
        Vector3 yAxis = new Vector3( 0.0F, 1.0F, 0.0F );

        //Compute bezel facet sides.
        Plane b = new Plane( lp[BT], lp[BB], new Vector3(lp[BT].x,lp[BT].y, -1.0F));
        float sScale = tableh*a1.x + (0.5F-tableh*a1.x)*star;
        lp[BR] = b.Intersect( new Vector3( sScale*a1.x, 0, -sScale*a1.y ),
                                  yAxis );

        //Compute pavilion facet sides.
        Plane p = new Plane( pb, lp[PT], new Vector3(pb.x,pb.y,-1.0F));
        float pScale = 0.5F*(1.0F-pavilion);
        lp[PR] = p.Intersect( new Vector3( pScale*a1.x, 0, -pScale*a1.y ),
                                  yAxis );

        //Compute girdle.
        float girdleRadius = 0.5F / (float)Math.cos( thetad / (girdleEdges*4) );
        Plane u0 = new Plane( lp[BR], lp[BB], lp[U0R] );
        Plane u1 = new Plane( lp[BR], lp[U0R], u1r );
        Plane l0 = new Plane( lp[PR], lp[PT], lp[L0R] );
        Plane l1 = new Plane( lp[PR], lp[L0R], l1r );
        for( int i = 0; i < girdleEdges; ++i )
        {
            Vector3 t0 = new Vector3( girdleRadius*(float)Math.cos(thetad/16.0*(2*i+1)),
                                     0.0F,
                                     -girdleRadius*(float)Math.sin(thetad/16.0*(2*i+1)));
            Vector3 t1 = new Vector3( girdleRadius*(float)Math.cos(thetad/16.0*(8+2*i+1)),
                                     0.0F,
                                     -girdleRadius*(float)Math.sin(thetad/16.0*(8+2*i+1)));
            lp[U0Iref+i] = u0.Intersect( t0, yAxis );
            lp[U1Iref+i] = u1.Intersect( t1, yAxis );
            lp[L0Iref+i] = l0.Intersect( t0, yAxis );
            lp[L1Iref+i] = l1.Intersect( t1, yAxis );

            //Sync
            lp[L0Iref+i].x = lp[U0Iref+i].x;
            lp[L0Iref+i].z = lp[U0Iref+i].z;
            lp[L1Iref+i].x = lp[U1Iref+i].x;
            lp[L1Iref+i].z = lp[U1Iref+i].z;
        }

        //Dupliate points
        for( int i = 0; i < 8; ++i )
        {
            Matrix3 r = Matrix3.MakeYRotation( thetad * i );
            for( int j = 0; j < ipoints; ++j )
            {
                points[i*ipoints+j] = new Vector3();
                Matrix3.Multiply( points[i*ipoints+j], r, lp[j] );
            }
        }
        if( !hasCulet )
        {
            points[8*ipoints] = new Vector3( pb );
        }

        //Place edges
        final int BT_BTp = 0;
        final int BT_BL = 1;
        final int BL_BB = 2;
        final int BB_BR = 3;
        final int BR_BT = 4;
        final int BR_U0R = 5;
        final int PT_PL = 6;
        final int PB_PR = 7;
        final int PR_PT = 8;
        final int PR_L0R = 9;
        final int U0Ieref = 10;
        final int U1Ieref = U0Ieref + girdleEdges+1;
        final int L0Ieref = U1Ieref + girdleEdges+1;
        final int L1Ieref = L0Ieref + girdleEdges+1;
        final int U0I_L0Iref = L1Ieref + girdleEdges+1;
        final int U1I_L1Iref = U0I_L0Iref + girdleEdges;
        final int C_Cp = (hasCulet ? iedges-1 : -1);

        Edge[] le = new Edge[iedges];
        
        le[BT_BTp] = new Edge( BT, BT+ipoints );
        le[BT_BL] = new Edge( BT, BR-ipoints );
        le[BL_BB] = new Edge( BR-ipoints, BB );
        le[BB_BR] = new Edge( BB, BR );
        le[BR_BT] = new Edge( BR, BT );
        le[BR_U0R] = new Edge( BR, U0R );
        le[PT_PL] = new Edge( PT, PR-ipoints );
        le[PB_PR] = new Edge( CL, PR );
        le[PR_PT] = new Edge( PR, PT );
        le[PR_L0R] = new Edge( PR, L0R );
        le[U0Ieref+0] = new Edge( BB, U0Iref+0 );
        le[U1Ieref+0] = new Edge( U0R, U1Iref+0 );
        le[L0Ieref+0] = new Edge( PT, L0Iref+0 );
        le[L1Ieref+0] = new Edge( L0R, L1Iref+0 );
        for( int i = 0; i < girdleEdges-1; ++i )
        {
            le[U0Ieref+1+i] = new Edge( U0Iref+i, U0Iref+i+1 );
            le[U1Ieref+1+i] = new Edge( U1Iref+i, U1Iref+i+1 );
            le[L0Ieref+1+i] = new Edge( L0Iref+i, L0Iref+i+1 );
            le[L1Ieref+1+i] = new Edge( L1Iref+i, L1Iref+i+1 );
        }
        le[U0Ieref+girdleEdges] = new Edge( U0Iref+girdleEdges-1, U0R );
        le[U1Ieref+girdleEdges] = new Edge( U1Iref+girdleEdges-1, 
                                                 BB+ipoints );
        le[L0Ieref+girdleEdges] = new Edge( L0Iref+girdleEdges-1, L0R );
        le[L1Ieref+girdleEdges] = new Edge( L1Iref+girdleEdges-1, 
                                                PT+ipoints );
        for( int i = 0; i < girdleEdges; ++i )
        {
            le[U0I_L0Iref+i] = new Edge( U0Iref+i, L0Iref+i );
            le[U1I_L1Iref+i] = new Edge( U1Iref+i, L1Iref+i );
        }

        //Duplicate edges
        for( int i = 0; i < 8; ++i )
        {
            for( int j = 0; j < iedges; ++j )
            {
                edges[i*iedges+j] = new Edge( (le[j].a+(8+i)*ipoints)%(8*ipoints),
                                              (le[j].b+(8+i)*ipoints)%(8*ipoints) );
            }
            if( !hasCulet )
            {
                edges[i*iedges+PB_PR].a = npoints-1;
            }
        }

        //Place faces
        final int B = 0;
        final int S = 1;
        final int U0 = 2;
        final int U1 = 3;
        final int P = 4;
        final int L0 = 5;
        final int L1 = 6;
        final int G0ref = 7;
        final int G1ref = G0ref + girdleEdges;
        final int C = (hasCulet ? ifaces-1 : -1);

        Face[] lf = new Face[ifaces];

        lf[B] = new Face( BT_BL, false, BL_BB, false, BB_BR, BR_BT );
        lf[S] = new Face( BR_BT, true, BT_BL+iedges, true, BT_BTp );
        lf[P] = new Face( PR_PT, false, PT_PL, false, PB_PR-iedges, PB_PR );
        int[] fu0 = new int[girdleEdges+1];
        int[] fu1 = new int[girdleEdges+1];
        int[] fl0 = new int[girdleEdges+1];
        int[] fl1 = new int[girdleEdges+1];
        for( int i = 0; i < girdleEdges+1; ++i )
        {
            fu0[i] = U0Ieref+i;
            fu1[i] = U1Ieref+i;
            fl0[girdleEdges-i] = L0Ieref+i;
            fl1[girdleEdges-i] = L1Ieref+i;
        }
        lf[U0] = new Face( BR_U0R, true, BB_BR, true, fu0 );
        lf[U1] = new Face( BL_BB+iedges, true, BR_U0R, false, fu1 );
        lf[L0] = new Face( PR_PT, true, PR_L0R, false, fl0 );
        lf[L1] = new Face( PR_L0R, true, PT_PL+iedges, true, fl1 );
        lf[G0ref+0] = new Face( L0Ieref+0, false, U0I_L0Iref+0, true,
                                U0Ieref+0, U1Ieref+girdleEdges - iedges,
                                U1I_L1Iref+girdleEdges-1 - iedges,
                                L1Ieref+girdleEdges - iedges );
        lf[G1ref+0] = new Face( L1Ieref+0, false, U1I_L1Iref+0, true,
                                U1Ieref+0, U0Ieref+girdleEdges,
                                U0I_L0Iref+girdleEdges-1,
                                L0Ieref+girdleEdges );
        for( int i = 1; i < girdleEdges; ++i )
        {
            lf[G0ref+i] = new Face( L0Ieref+i, false, U0I_L0Iref+i, true,
                                    U0Ieref+i, U0I_L0Iref+i-1 );
            lf[G1ref+i] = new Face( L1Ieref+i, false, U1I_L1Iref+i, true,
                                    U1Ieref+i, U1I_L1Iref+i-1 );
        }

        //Duplicate faces
        for( int i = 0; i < 8; ++i )
        {
            for( int j = 0; j < ifaces; ++j )
            {
                int[] r = new int[lf[j].edges.length];
                for( int k = 0; k < r.length; ++k )
                {
                    r[k] = (lf[j].edges[k] + (8+i)*iedges)%(8*iedges);
                }
                faces[i*ifaces+j] = new Face( r, lf[j].reverse0,
                                              lf[j].reverse1 );
            }
        }
        int[] table = new int[8];
        for( int i = 0; i < 8; ++i )
        {
            table[i] = BT_BTp + i*iedges;
        }
        faces[nfaces-1] = new Face( table, false, false );

        Parameters retval = new Parameters();
        retval.points = points;
        retval.edges = edges;
        retval.faces = faces;

        return retval;
    }

    static private class Plane
    {
        Vector3 ref;
        Vector3 norm;

        Plane( Vector3 p1, Vector3 p2, Vector3 p3 )
        {
            ref = new Vector3( p1 );

            Vector3 l0 = new Vector3();
            Vector3.Minus( l0, p2, p1 );
            Vector3 l1 = new Vector3();
            Vector3.Minus( l1, p3, p1 );
            norm = new Vector3();
            Vector3.Cross( norm, l0, l1 );
            Vector3.Normalize( norm );
        }

        Vector3 Intersect( Vector3 iref, Vector3 idir )
        {
            Vector3 i = new Vector3();
            Vector3 rel = new Vector3();
            Vector3.Minus( rel, iref, ref );
            Vector3.IntersectPlane( i, norm, rel, idir );
            Vector3.Add( i, i, ref );
            return i;
        }
    }

    static RoundCut MakeTest()
    {
        //1.26 carat
        Measurements m = new Measurements();
        m.minWidth = new Float( 0.0690F );
        m.maxWidth = new Float( 0.0694F );
        m.height = new Float( 0.0426F );
        m.tablePct = new Float( 0.57F );
        //m.heightPct = new Float( 0.615F );
        //m.crownHeightPct; = new Float( 0.15F );
        //m.pavilionHeightPct; = new Float( 0.445F );
        m.crownAngleDeg = new Float( 34.5F );
        m.pavilionAngleDeg = new Float( 40.8F );
        //m.girdlePct; = new Float( 0.035F );
        m.starPct = new Float( 0.55F );
        m.pavilionPct = new Float( 0.75F );
        m.girdleEdges = new Integer( 4 );
        m.hasCulet = false;
        m.Complete();
        return new RoundCut( m );
    }
}
