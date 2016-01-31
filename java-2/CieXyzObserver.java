import java.io.File;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.io.StringReader;
import java.util.Scanner;

public class CieXyzObserver {
    static final int MinWavelength = 380;
    static final int MaxWavelength = 780;
    static final int Increment = 5;
    
    static private class Datum
    {
        float x;
        float y;
        float z;
	
        Datum()
        {
            x = 0;
            y = 0;
            z = 0;
        }
	
        Datum( Datum that )
        {
            x = that.x;
            y = that.y;
            z = that.z;
        }
	
        public Datum Add( Datum that )
        {
            Datum retval = new Datum();
            retval.x = x + that.x;
            retval.y = y + that.y;
            retval.z = z + that.z;
            
            return retval;
        }	
	
        public Datum Multiply( float s )
        {
            Datum retval = new Datum();
            retval.x = x * s;
            retval.y = y * s;
            retval.z = z * s;
            
            return retval;
        }
    }
    
    Datum[] mData;
    
    public CieXyzObserver( String filename )
	throws java.io.FileNotFoundException, java.io.IOException
    {
        mData = new Datum[(MaxWavelength-MinWavelength)/Increment + 1];
	
        File f = new File( filename );
        FileReader fr = new FileReader( f );
        LineNumberReader in = new LineNumberReader( fr );
	
        while(true)
        {
            
            String line = in.readLine();
            if( null != line )
            {
                StringReader sline = new StringReader( line );
                Scanner s = new Scanner( sline );
                s.useDelimiter(",");
                int index = GetIndex( s.nextInt() );
                if( index+1 == in.getLineNumber() )
                {
                    mData[index] = new Datum();
                    mData[index].x = s.nextFloat();
                    mData[index].y = s.nextFloat();
                    mData[index].z = s.nextFloat();
                }
                else
                {
                    throw new java.io.IOException( "Bad Format" );
                }
            }
            else
            {
                break;
            }
        }
	
        for(int i = 0; i < mData.length; ++i )
        {
            if( null == mData[i] )
            {
                throw new java.io.IOException( "Missing data" );
            }
        }
    }
    
    public Datum GetDatum( float wavelength )
    {
        int index = GetIndex( wavelength );
        if( index != -1 )
        {
            int index2;
            if( index < mData.length - 1 )
            {
                index2 = index + 1;
            }
            else
            {
                index2 = index;
            }
			
            Datum retval = new Datum();
            float p = (wavelength - (MinWavelength+index*Increment))/Increment;
            retval.x = mData[index].x + p * (mData[index2].x - mData[index].x);
            retval.y = mData[index].y + p * (mData[index2].y - mData[index].y);
            retval.z = mData[index].z + p * (mData[index2].z - mData[index].z);
            return retval;
        }
        else
        {
            return null;
        }
    }

    public CieXyz Integrate( float w0, float w1 )
    {
        CieXyzObserver.Datum datum;
        if( w0 < MaxWavelength && w1 > MinWavelength )
        {
            if( w0 < MinWavelength ) w0 = MinWavelength;
            if( w1 > MaxWavelength ) w1 = MaxWavelength;
		
            int i0 = GetIndex( w0 );
            int i1 = GetIndex( w1 );
			
            if( i0 == i1 )
            {
                datum = new CieXyzObserver.Datum( mData[i0] ).Multiply( (w1-w0)/(float)Increment );
            }
            else
            {
                datum = new CieXyzObserver.Datum();
                datum = datum.Add( mData[i0].Multiply( ((MinWavelength+(i0+1)*Increment) - w0)/(float)Increment ) );
                for( int i = i0+1; i < i1; ++i )
                {
                    datum = datum.Add( mData[i] );
                }
                datum = datum.Add( mData[i1].Multiply( (w1 - (MinWavelength+(i1)*Increment))/(float)Increment ) );
            }
        }
        else
        {
            datum = new CieXyzObserver.Datum();
        }
        datum = datum.Multiply( 1.0F / 21.3715F );
        CieXyz xyz = new CieXyz();
        xyz.X = Bound( datum.x );
        xyz.Y = Bound( datum.y );
        xyz.Z = Bound( datum.z );
		
        return xyz;
    }
	
    private int GetIndex( float wavelength )
    {
        int iwave = (int)wavelength;
        if( iwave >= MinWavelength && iwave <= MaxWavelength )
        {
            return (iwave-MinWavelength)/Increment;
        }
        else
        {
            return -1;
        }
    }
	
    private float Bound( float s )
    {
        if( s > 1.0 )
        {
            return 1.0F;
        }
        else if( s >= 0.0 )
        {
            return s;
        }
        else
        {
            return 0.0F;
        }
    }
}
