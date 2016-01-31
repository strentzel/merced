import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*;

class ToPng
{
    static public void main(String[] args) throws IOException
    {
        Scanner in = new Scanner(System.in);
        int width = in.nextInt();
        int height = in.nextInt();
        int[] buffer = new int[width*height];
        for (int i = 0; i < buffer.length; ++i) {
            int a = in.nextInt();
            int r = in.nextInt();
            int g = in.nextInt();
            int b = in.nextInt();
            int argb = (a << 24) | (r << 16) | (g << 8) | b;
            buffer[i] = argb;
        }

        BufferedImage image = new BufferedImage( width, height, BufferedImage.TYPE_INT_ARGB );
        image.setRGB( 0, 0, width, height, buffer, 0, width );

        ImageIO.write( image, args[0], System.out );
    }
}
