public interface IDrawable
{
    static public float LineWidth = 3.0F;
    public void Rasterize( float scale, Vector2 transform, 
                           RawBitmap bitmap, int rgb );
}
