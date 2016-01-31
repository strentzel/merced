public class Material
{
    public final float MinN;
    public final float MaxN;
    public final float Theta;
    public final float CosTheta;
    private final float C0;
    private final float C1;
    private final float C2;

    public Material( float a, float b, float c )
    {
        C0 = -b / (2.0F * a);
        C1 = 1.0F / a;
        C2 = (b / (2.0F * a)) * (b / (2.0F * a)) - c / a; 
        MinN = -C2 / C1;
        MaxN = (float) (a*380*380+b*380+c);
        Theta = (float) Math.asin( 1.0 / MinN );
        CosTheta = (float) Math.cos( Theta );
    }

    public Material( Material that )
    {
        C0 = that.C0;
        C1 = that.C1;
        C2 = that.C2;
        MinN = that.MinN;
        MaxN = that.MaxN;
        Theta = that.Theta;
        CosTheta = that.CosTheta;
    }

    public float EstimateWavelength( float n )
    {
        if( n > MinN )
        {
            return C0 - (float)Math.sqrt(C1 * n + C2);
        }
        else
        {
            return C0;
        }
    }

    public float EstimateReflectionWeight( float thetaIncident )
    {
        float thetaTransmittance = (float)Math.asin( 
            Math.sin( thetaIncident ) * MinN );
        float cost = (float)Math.cos( thetaTransmittance );
        float cosi = (float)Math.cos( thetaIncident );
        float Rs = (float)Math.pow( (MinN*cosi - cost)/
                                    (MinN*cosi + cost), 2.0 );
        float Rp = (float)Math.pow( (MinN*cost - cosi)/
                                    (MinN*cost + cosi), 2.0 );
        return (Rs+Rp)/2.0F;
    }

    static public Material GetDiamond()
    {
        return new Material( 4.98E-7F, -7.28E-4F, 2.67F );
    }
}
