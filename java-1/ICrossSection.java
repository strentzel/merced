public interface ICrossSection
{
    public int GetCrossSectionCnt();
    public Matrix3 GetCrossSectionRot( int index );
    public Vector3[] GetCrossSection( int index );
}
