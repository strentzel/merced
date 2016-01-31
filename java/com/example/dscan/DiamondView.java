package com.example.dscan;

import java.util.ListIterator;

import android.content.Context;
import android.util.AttributeSet;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Path;
import android.view.MotionEvent;
import android.view.View;

public class DiamondView extends View {

	RoundCut mDiamond;
	Material mInternal;
	Material mExternal;
	Observer mObserver;
	ConvexPolygon mCrossSection;

	AndroidCanvas mCanvas;
	int[] mBitmap;
	
	private boolean mStateDown;
	private float mStateX0;
	private float mStateX1;
	private float mStateY0;
	private float mStateY1;
	
	Vector2 mVec1;
	Vector2 mVec2;
	Vector2 mVec3;
	Vector2 mVec4;
	Matrix3 mMat1;
	Matrix3 mMat2;
	
	Path mPath;

	public DiamondView(Context context) {
		this(context, null, 0);
	}

	public DiamondView(Context context, AttributeSet attrs) {
		this(context, attrs, 0);
	}

	public DiamondView(Context context, AttributeSet attrs, int defStyle) {
		super(context, attrs, defStyle);

		mCanvas = new AndroidCanvas();
		
		mStateX0 = 0;
		mStateY0 = 0;
//		mStateX1 = 78;
//		mStateY1 = 149;
		mStateX1 = 24;
		mStateY1 = 172;
		
                RoundCut.Measurements m = new RoundCut.Measurements();
                m.tablePct = 0.57F;
                m.crownHeightPct = 0.15F;
                m.pavilionHeightPct = 0.445F;
                m.girdlePct = 0.035F;
                m.starPct = 0.55F;
                m.pavilionPct = 0.75F;
                m.culetPct = 0.0F;
                m.girdleEdges = 4;

		mDiamond = new RoundCut(m);
                mInternal = new Material( "diamond" );
                mExternal = new Material( "air" );
                mObserver = new Observer( "default" );
                mCrossSection = mDiamond.GetCrossSection(0);
		
		mVec1 = new Vector2();
		mVec2 = new Vector2();
		mVec3 = new Vector2();
		mVec4 = new Vector2();
		mMat1 = new Matrix3();
		mMat2 = new Matrix3();
		
		mPath = new Path();
	}

	@Override
	public boolean onTouchEvent(MotionEvent event) {
		switch (event.getActionMasked()) {
		case MotionEvent.ACTION_DOWN:
			mStateDown = true;
			mStateX0 = mStateX1 = event.getX();
			mStateY0 = mStateY1 = event.getY();
			return true;
		case MotionEvent.ACTION_MOVE:
			mStateX1 = event.getX();
			mStateY1 = event.getY();
			invalidate();
			return true;
		case MotionEvent.ACTION_CANCEL:
			mStateDown = false;
			invalidate();
			return true;
		case MotionEvent.ACTION_UP:
			mStateDown = false;
			mStateX1 = event.getX();
			mStateY1 = event.getY();
			invalidate();
			return true;
		default:
			return false;
		}
	}
	
	@Override
	public void onDraw(Canvas canvasImpl) {
            Native.EnableExceptions();

            mCanvas.set(canvasImpl);
		
            int width = getMeasuredWidth();
            int height = getMeasuredHeight();
            
            Matrix3 displayTransform = mMat1;
            float scale = (float)Math.min(width, height) * 2.0F / 3.0F; 
            displayTransform.d00 = scale;
            displayTransform.d02 = width/2.0F;
            displayTransform.d11 = -scale;
            displayTransform.d12 = height/2.0F;
            displayTransform.d22 = 1.0F;
            
            Edge[] edges = mDiamond.GetVisibleEdges();
            for (int i = 0; i < edges.length; ++i) {
                mCanvas.drawLine( 
                                 displayTransform.homogenousTransform(edges[i].start),
                                 displayTransform.homogenousTransform(edges[i].end),
                                 0xFF007F7F, 3);
            }
            Vector2 manual = mVec3;
            manual.x = mStateX1;
            manual.y = mStateY1;
            Matrix3 idt = mMat2;
            Matrix3.homogenousInverse(idt, displayTransform);
            manual = idt.homogenousTransform(manual);
            
            LightPath path = new LightPath( new Vector2( 0.0F, 0.5F ),
                                            manual, mExternal, mInternal,
                                            mObserver, mCrossSection,
                                            1E-4F, 1E-4F );
            ListIterator<LightPath.Step> sit = path.GetSteps();
            while( sit.hasNext() ){
                LightPath.Step step = sit.next();
                
                for( int i = 0; i < step.points.length; ++i ){
                    Vector2 point = displayTransform.homogenousTransform( step.points[i] );
                    if (i == 0 ) {
                        mPath.moveTo(point.x, point.y);
                    } else {
                        mPath.lineTo(point.x, point.y);
                    }
                }
                mPath.close();
                mCanvas.drawPath(mPath, Color.argb((int)(64.0*step.weight), 255, 255, 255));
            }
		
            for (int i = 0; i < width*height; ++i) {
                mBitmap[i] = 0;
            }

            ListIterator<Prism> prisms = path.GetPrisms();
            while( prisms.hasNext() ){
                Prism prism = prisms.next();
                prism.Rasterize( displayTransform, 1.0F, mBitmap, 0, width, width, height );
            }
	
            mCanvas.drawBitmap(mBitmap, width, height);
		
            Native.DisableExceptions();
        }

	@Override
	protected void onMeasure(int widthMeasureSpec, int heightMeasureSpec) {
		setMeasuredDimension(widthMeasureSpec, heightMeasureSpec);
	}

	@Override
	protected void onSizeChanged(int width, int height, int oldWidth, int oldHeight) {
		mBitmap = new int[width*height];
	}
}
