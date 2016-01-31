public class Material
{
    static private final float DiamondMinN = 4834000.0F / 2008000.0F;
    static private final float DiamondMaxN = (float)(((731.3-380)*(731.3-380)+4834000.0)/2008000.0);

    public Material( float mn, float mx )
    {
        MinN = mn;
        MaxN = mx;
        Theta = (float) Math.asin( 1 / MinN );
        CosTheta = (float)Math.cos( Theta );
    }

    public Material( Material that )
    {
        MinN = that.MinN;
        MaxN = that.MaxN;
        Theta = that.Theta;
        CosTheta = that.CosTheta;
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
        return new Material( DiamondMinN, DiamondMaxN );
    }

    public final float MinN;
    public final float MaxN;
    public final float Theta;
    public final float CosTheta;
}
