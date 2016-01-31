package com.example.dscan;

import android.annotation.TargetApi;
import android.graphics.Paint;
import android.graphics.Path;
import android.graphics.PorterDuff;
import android.graphics.PorterDuffXfermode;

@TargetApi(11)
public class AndroidCanvas {
	private android.graphics.Canvas mAndroidCanvas;
	private Paint mPixelPaint;
	private Paint mLinePaint;
	private Paint mShapePaint;
	
	public AndroidCanvas() {
		mAndroidCanvas = null;
		mPixelPaint = new Paint();
		mPixelPaint.setXfermode(new PorterDuffXfermode(PorterDuff.Mode.ADD));
		mLinePaint = new Paint(Paint.ANTI_ALIAS_FLAG);
		mLinePaint.setStyle(Paint.Style.STROKE);
		mLinePaint.setXfermode(new PorterDuffXfermode(PorterDuff.Mode.ADD));
		mShapePaint = new Paint(Paint.ANTI_ALIAS_FLAG);
		mShapePaint.setStyle(Paint.Style.FILL);
		mShapePaint.setXfermode(new PorterDuffXfermode(PorterDuff.Mode.ADD));
	}

	public void set( android.graphics.Canvas canvas ) {
		mAndroidCanvas = canvas;
	}
	
	public void drawBitmap(int[] bitmap, int width, int height) {
		mAndroidCanvas.drawBitmap(bitmap, 0, width, 0, 0, width, height, true, mPixelPaint);
	}
	
	public void drawLine(Vector2 start, Vector2 end, int argb, float width) {
		mLinePaint.setStrokeWidth(width);
		mLinePaint.setColor(argb);
		mAndroidCanvas.drawLine(start.x, start.y, end.x, end.y, mLinePaint);
	}

	public void drawPath(Path path, int argb) {
		mShapePaint.setColor(argb);
		mAndroidCanvas.drawPath(path, mShapePaint);
	}
}
