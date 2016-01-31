
public class CieXyz {
	public float X;
	public float Y;
	public float Z;
	
	public CieXyz AddMult( CieXyz that, float s )
	{
		CieXyz retval = new CieXyz();
		retval.X = X + that.X * s;
		retval.Y = Y + that.Y * s;
		retval.Z = Z + that.Z * s;
		return retval;
	}
	
	public int toRgb( int bg )
	{
		float rf =  2.364614F*X - 0.896541F*Y - 0.468073F*Z;
		float gf = -0.515166F*X + 1.426408F*Y + 0.088758F*Z;
		float bf =  0.005204F*X - 0.014408F*Y + 1.009204F*Z;
		float mx = rf;
		if( gf > mx ) mx = gf;
		if( bf > mx ) mx = bf; 
		float Yp = (Y<1.0F?Y:1.0F);
		int r = cap( 255 * Yp * rf / mx + (bg/(256*256)) * (1-Y) );
		int g = cap( 255 * Yp * gf / mx + ((bg/256)%256) * (1-Y) );
		int b = cap( 255 * Yp * bf / mx + (bg%256) * (1-Y) );
		return r * 256 * 256 + g * 256 + b;
	}
	
	private int cap( float f )
	{
		if( f < 0 ) return 0;
		else if( f > 254.5F ) return 255;
		else return (int)(f+0.5F);
	}
}
