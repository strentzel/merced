import gameoxide.GOImageCanvas;

public class RawBitmap
{
    int width;
    int height;
    int[] data;

    public RawBitmap( int w, int h )
    {
        width = w;
        height = h;
        data = new int[width*height];
    }

    public void Add( int x, int y, int r, int g, int b )
    {
        if( x >= 0 && x < width && y >= 0 && y < height )
        {
            int curr = data[y*width + x];
            int rp = ((curr >> 16) % 256) + r;
            if( rp > 255 ) rp = 255;
            int gp = ((curr >> 8) % 256) + g;
            if( gp > 255 ) gp = 255;
            int bp = (curr % 256) + b;
            if( bp > 255 ) bp = 255;
            data[y*width + x] = (rp << 16) + (gp << 8) + bp;
        }
    }

    public void DrawTriangle( Vector2 a, Vector2 b, Vector2 c, int rgb )
    {
        Vector2 e1 = new Vector2();
        Vector2 e2 = new Vector2();
        Vector2 e3 = new Vector2();
        Vector2.Minus( e1, b, a );
        Vector2.Minus( e2, c, b );
        Vector2.Minus( e3, a, c );
        Vector2 n1 = new Vector2();
        Vector2 n2 = new Vector2();
        Vector2 n3 = new Vector2();
        Vector2.NormFromRel( n1, e1, e2 );
        Vector2.NormFromRel( n2, e2, e3 );
        Vector2.NormFromRel( n3, e3, e1 );
        if( n1.Dot( n2 ) < 0 || n2.Dot( n3 ) < 0 || n3.Dot( n1 ) < 0 )
        {
            System.out.println( n1 + " " + n2 + " " + n3 );
        }

        Vector2 p = new Vector2();
        Vector2 tmp = new Vector2();
        for( int row = 0; row < height; ++row )
        {
            p.y = row;
            for( int col = 0; col < width; ++col )
            {
                p.x = col;
                Vector2.Minus( tmp, p, a );
                if( tmp.Dot( n1 ) >= 0 )
                {
                    Vector2.Minus( tmp, p, b );
                    if( tmp.Dot( n2 ) >= 0 )
                    {
                        Vector2.Minus( tmp, p, c );
                        if( tmp.Dot( n3 ) >= 0 )
                        {
                            int x = data[row*width+col];
                            data[row*width+col] =
                                (Math.min(255,
                                          ((x>>16)%256)+
                                          ((rgb>>16)%256))<<16)+
                                (Math.min(255,
                                          ((x>>8)%256)+
                                          ((rgb>>8)%256))<<8)+
                                (Math.min(255,
                                          ((x>>0)%256)+
                                          ((rgb>>0)%256))<<0);
                        }
                    }
                }
            }
        }
    }

    public void Draw( GOImageCanvas canvas )
    {
        int index = 0;
        for( int row = 0; row < height; ++row )
        {
            for( int col = 0; col < width; ++col )
            {
                canvas.moveTo( col, row );
                canvas.setColor( data[row*width+col] );
                canvas.drawTo( col, row+1 );
                ++index;
            }
        }
    }
}