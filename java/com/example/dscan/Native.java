package com.example.dscan;

public class Native {
	native static public void EnableExceptions();
	native static public void DisableExceptions();
        native static public long RoundCutCreate(float[] values, int girdleEdges);
	native static public int RoundCutGetRotations();
	native static public float[] RoundCutGetRotation(int r);
	native static public long RoundCutGetCrossSection(long poly, int index);
	native static public long MaterialCreate(String name);
	native static public void MaterialDestroy(long material);
	native static public long ObserverCreate(String name);
	native static public void ObserverDestroy(long observer);
	native static public void ConvexPolyhedronDestroy(long poly);
	native static public void ConvexPolyhedronTransform(long poly, float d00, float d01, float d02, float d10, float d11, float d12, float d20, float d21, float d22);
	native static public float ConvexPolyhedronVolume(long poly);
	native static public float[] ConvexPolyhedronVisibleEdges(long poly);

        native static public void ConvexPolygonDestroy(long poly);
        native static public boolean ConvexPolygonTrace(long path, long external, long internal, long crossSection, float minWeight, float minSurface);
        native static public long LightPathCreate(float sx, float sy, float tx, float ty);
	native static public void LightPathDestroy(long path);
        native static public int LightPathElements(long path);
        native static public float[] LightPathMarshall(long path, int first, int last);

        native static public long PrismCreate(float sx, float sy, float tx, float ty, float ix, float iy, long materialIncident, long materialTransmittent, long observer, float d00, float d01, float d02, float d10, float d11, float d12, float weight);
        native static public void PrismDestroy(long prism);
        native static public void PrismRasterize(long prism, int[] bitmap, int offset, int stride, int width, int height);
}
