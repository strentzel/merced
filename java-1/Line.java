class Line
{
    Triangle[] polys;

    public Line( Vector2 a, Vector2 b, float width )
    {
        Vector2 norm = new Vector2( b.y - a.y, a.x - b.x );
        Vector2.Normalize( norm );

        float n2 = norm.Norm2();
        if( !(n2 > 0.9F && n2 < 1.1F) )
        {
            norm.x = 0;
            norm.y = 0;
        }
        
        Vector2[] points = new Vector2[4];
        points[0] = new Vector2( norm );
        Vector2.Scale( points[0], points[0], width/2.0F );
        Vector2.Add( points[0], points[0], a );
        points[1] = new Vector2( norm );
        Vector2.Scale( points[1], points[1], -width/2.0F );
        Vector2.Add( points[1], points[1], a );
        points[2] = new Vector2( norm );
        Vector2.Scale( points[2], points[2], width/2.0F );
        Vector2.Add( points[2], points[2], b );
        points[3] = new Vector2( norm );
        Vector2.Scale( points[3], points[3], -width/2.0F );
        Vector2.Add( points[3], points[3], b );

        polys = new Triangle[2];
        polys[0] = new Triangle( points[0], points[1], points[2] );
        polys[1] = new Triangle( points[1], points[2], points[3] );
    }

    public void Rasterize( RawBitmap bitmap, int rgb )
    {
        polys[0].Rasterize( bitmap, rgb );
        polys[1].Rasterize( bitmap, rgb );
    }
}
