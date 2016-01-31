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

    public void Add( int x, int y, int rgb )
    {
        Add( x, y, (rgb >> 16)%256, (rgb >> 8)%256, rgb%256 );
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
