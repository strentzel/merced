import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

class LightBeam
{
    static public class CrossSection
    {
        public Vector2 start;
        public Vector2 end;
        public float weight;

        public CrossSection scale( float s, float t, Vector2 work )
        {
            CrossSection retval = new CrossSection();
            retval.start = new Vector2();
            retval.end = new Vector2();
            retval.weight = weight;
            Vector2.Minus( work, end, start );
            Vector2.Scale( retval.start, work, s );
            Vector2.Add( retval.start, start, retval.start );
            Vector2.Scale( retval.end, work, t );
            Vector2.Add( retval.end, start, retval.end );
            return retval;
        }            

        public String toString()
        {
            return start + "to" + end + "(" + weight + ")";
        }
    }

    List<CrossSection> scalePath( float s, float t, Vector2 work )
    {
        List<CrossSection> retval = new LinkedList<CrossSection>();
        ListIterator<CrossSection> it = path.listIterator();
        while( it.hasNext() )
        {
            retval.add( it.next().scale( s, t, work ) );
        }
        retval.add( src.scale( s, t, work ) );
        return retval;
    }

    public void Draw( RawBitmap bitmap )
    {
        if( !path.isEmpty() )
        {
            ListIterator<CrossSection> it = path.listIterator();
            CrossSection prev = it.next();
            while( it.hasNext() )
            {
                int rgb = (((int)(255*prev.weight))<<16) +
                    (((int)(255*prev.weight))<<8) +
                    ((int)(255*prev.weight));
                CrossSection curr = it.next();
//                bitmap.DrawTriangle( prev.start, prev.end, curr.start, rgb );
//                bitmap.DrawTriangle( curr.start, curr.end, prev.end, rgb );
                Triangle t1 = new Triangle( prev.start, prev.end, curr.start );
                t1.Rasterize( bitmap, rgb );
                Triangle t2 = new Triangle( curr.start, curr.end, prev.end );
                t2.Rasterize( bitmap, rgb );
                prev = curr;
            }
            int rgb = (((int)(255*prev.weight))<<16) +
                (((int)(255*prev.weight))<<8) +
                ((int)(255*prev.weight));
//            bitmap.DrawTriangle( prev.start, prev.end, src.start, rgb ); 
//            bitmap.DrawTriangle( src.start, src.end, prev.end, rgb ); 
            Triangle t1 = new Triangle( prev.start, prev.end, src.start );
            t1.Rasterize( bitmap, rgb );
            Triangle t2 = new Triangle( src.start, src.end, prev.end );
            t2.Rasterize( bitmap, rgb );
        }

        CieXyzObserver observer = null;
        try
        {
            observer = new CieXyzObserver( "data/cie.csv" );
        }
        catch( java.io.IOException ex ){System.exit(0);}
        Prism P = new Prism( observer, src.start, src.end, dir, weight );
        Vector2 tmp1 = new Vector2();
        for( int row = 0; row < 480; ++row )
        {
            tmp1.y = row;
            for( int col = 0; col < 640; ++col )
            {
                tmp1.x = col;
                int rgb = P.Compute( tmp1 );
                if( rgb > 0 )
                {
                    bitmap.Add( col, row, ((rgb>>16)%256),
                                ((rgb>>8)%256), (rgb%256) );
                }
            }
        }        
    }

    List<CrossSection> path;
    public CrossSection src;
    public Vector2 dir;
    public float weight;
    public boolean inside;

    public String toString()
    {
        String retval = "LB: " + (int)(weight*100) + "%";
        ListIterator<CrossSection> it = path.listIterator();
        while( it.hasNext() )
        {
            CrossSection cs = it.next();
            retval += " " + cs;
        }
        retval += " " + src;
        retval += " + " + dir;
        return retval;
    }
}