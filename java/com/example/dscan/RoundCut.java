package com.example.dscan;

public class RoundCut extends ConvexPolyhedron
{
    static public class Measurements
    {
        public float tablePct;
        public float crownHeightPct;
        public float pavilionHeightPct;
        public float girdlePct;
        public float starPct;
        public float pavilionPct;
        public float culetPct;
        public int girdleEdges;
    }

    static private float[] MakeValues( Measurements m )
    {
        float[] values = new float[7];
        values[0] = m.tablePct;
        values[1] = m.crownHeightPct;
        values[2] = m.pavilionHeightPct;
        values[3] = m.girdlePct;
        values[4] = m.starPct;
        values[5] = m.pavilionPct;
        values[6] = m.culetPct;
        return values;
    }

    public RoundCut( Measurements m )
    {
        super( Native.RoundCutCreate( MakeValues(m), m.girdleEdges ) );
    }

    public int GetRotations()
    {
        return Native.RoundCutGetRotations();
    }

    public Matrix3 GetRotation( int r )
    {
        float[] values = Native.RoundCutGetRotation(r);
        Matrix3 retval = new Matrix3();
        retval.d00 = values[0];
        retval.d01 = values[1];
        retval.d02 = values[2];
        retval.d10 = values[3];
        retval.d11 = values[4];
        retval.d12 = values[5];
        retval.d20 = values[6];
        retval.d21 = values[7];
        retval.d22 = values[8];
        return retval;
    }

    public ConvexPolygon GetCrossSection( int r )
    {
        return new ConvexPolygon( Native.RoundCutGetCrossSection( mImpl, r ) );
    }
}
