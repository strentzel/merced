package com.example.dscan;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

public class LightPath
{
    private class Element
    {
        int type;
        float weight;
        Vector2[] points;        
    }

    public class Step
    {
        float weight;
        Vector2[] points;

        Step( Element e )
        {
            weight = e.weight;
            points = e.points;
        }
    }

    public LightPath( Vector2 start, Vector2 end, Material external,
                      Material internal, Observer observer, 
                      ConvexPolygon crossSection, float minWeight,
                      float minSurface )
    {
        long path = Native.LightPathCreate( start.x, start.y, end.x, end.y );
        for( int i = 0; i < 50; ++i ) {
            if( !Native.ConvexPolygonTrace( path, external.mImpl,
                                            internal.mImpl, crossSection.mImpl,
                                            minWeight, minSurface ) ){
                break;
            }
        }

        float[] values = Native.LightPathMarshall( path, 0, 
                                                   Native.LightPathElements(path));
        int index = 0;
        int elementsCnt = (int) values[index++];
        int[] cnts = new int[2];
        Element[] elements = new Element[elementsCnt];
        for (int i = 0; i < elementsCnt; ++i) {
            elements[i] = new Element();
            elements[i].type = (int) values[index++];
            elements[i].weight = values[index++];
            int cnt = (int) values[index++];
            elements[i].points = new Vector2[cnt];
            for (int j = 0; j < cnt; ++j) {
                elements[i].points[j] = new Vector2();
                elements[i].points[j].x = values[index++];
                elements[i].points[j].y = values[index++];
            }
            cnts[elements[i].type]++;
        }

        mSteps = new ArrayList<Step>( cnts[0] );
        mPrisms = new ArrayList<Prism>( cnts[1] );

        for( int i = 0; i < elements.length; ++i ){
            if( elements[i].type == 0 ){
                mSteps.add( new Step( elements[i] ) );
            } else {
                mPrisms.add( new Prism( external,
                                        internal,
                                        observer,
                                        elements[i].points[0], 
                                        elements[i].points[1],
                                        elements[i].points[2],
                                        elements[i].weight ) );
            }
        }

        Native.LightPathDestroy( path );
    }

    public ListIterator<Step> GetSteps()
    {
        return mSteps.listIterator();
    }

    public ListIterator<Prism> GetPrisms()
    {
        return mPrisms.listIterator();
    }

    private List<Step> mSteps;
    private List<Prism> mPrisms;
}
