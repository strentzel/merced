package com.example.dscan;

class Matrix3 {
	public float d00;
	public float d01;
	public float d02;
	public float d10;
	public float d11;
	public float d12;
	public float d20;
	public float d21;
	public float d22;

	public Matrix3() {
		d00 = d01 = d02 = d10 = d11 = d12 = d20 = d21 = d22 = 0.0F;
	}

	public Vector2 homogenousTransform(Vector2 point) {
            Vector2 result = new Vector2();
		assert d20 == 0.0 && d21 == 0.0 && d22 == 1.0;
		float rx = d00 * point.x + d01 * point.y + d02;
		float ry = d10 * point.x + d11 * point.y + d12;
		result.x = rx;
		result.y = ry;
                return result;
	}
	
	static public Matrix3 MakeIdentity() {
		Matrix3 I = new Matrix3();
		I.d00 = 1.0F;
		I.d11 = 1.0F;
		I.d22 = 1.0F;
		return I;
	}

	static public Matrix3 MakeXRotation(float theta) {
		Matrix3 r = new Matrix3();
		r.d00 = 1.0F;
		r.d10 = 0.0F;
		r.d20 = 0.0F;
		r.d01 = 0.0F;
		r.d11 = (float) NearCos(theta);
		r.d21 = (float) NearSin(theta);
		r.d02 = 0.0F;
		r.d12 = (float) -NearSin(theta);
		r.d22 = (float) NearCos(theta);

		return r;
	}

	static public Matrix3 MakeYRotation(float theta) {
		Matrix3 r = new Matrix3();
		r.d00 = (float) NearCos(theta);
		r.d10 = 0.0F;
		r.d20 = (float) -NearSin(theta);
		r.d01 = 0.0F;
		r.d11 = 1.0F;
		r.d21 = 0.0F;
		r.d02 = (float) NearSin(theta);
		r.d12 = 0.0F;
		r.d22 = (float) NearCos(theta);

		return r;
	}

	static public Matrix3 MakeZRotation(float theta) {
		Matrix3 r = new Matrix3();
		r.d00 = (float) NearCos(theta);
		r.d10 = (float) NearSin(theta);
		r.d20 = 0.0F;
		r.d01 = (float) -NearSin(theta);
		r.d11 = (float) NearCos(theta);
		r.d21 = 0.0F;
		r.d02 = 0.0F;
		r.d12 = 0.0F;
		r.d22 = 1.0F;

		return r;
	}

	static public void Multiply(Matrix3 r, Matrix3 m, float s) {
		r.d00 = m.d00 * s;
		r.d01 = m.d01 * s;
		r.d02 = m.d02 * s;
		r.d10 = m.d10 * s;
		r.d11 = m.d11 * s;
		r.d12 = m.d12 * s;
		r.d20 = m.d20 * s;
		r.d21 = m.d21 * s;
		r.d22 = m.d22 * s;
	}

	static public void homogenousInverse(Matrix3 r, Matrix3 h) {
		float idet = 1.0F / (h.d00*h.d11 - h.d01*h.d10);
		r.d00 = h.d11 * idet;
		r.d01 = -h.d01 * idet;
		r.d10 = -h.d10 * idet;
		r.d11 = h.d00 * idet;
		r.d02 = -(r.d00*h.d02 + r.d01*h.d12);
		r.d12 = -(r.d10*h.d02 + r.d11*h.d12);
		r.d20 = 0.0F;
		r.d21 = 0.0F;
		r.d22 = 1.0F;
	}
	
	static public void Multiply(Matrix3 r, Matrix3 a, Matrix3 b) {
		float r00 = a.d00 * b.d00 + a.d01 * b.d10 + a.d02 * b.d20;
		float r01 = a.d00 * b.d01 + a.d01 * b.d11 + a.d02 * b.d21;
		float r02 = a.d00 * b.d02 + a.d01 * b.d12 + a.d02 * b.d22;
		float r10 = a.d10 * b.d00 + a.d11 * b.d10 + a.d12 * b.d20;
		float r11 = a.d10 * b.d01 + a.d11 * b.d11 + a.d12 * b.d21;
		float r12 = a.d10 * b.d02 + a.d11 * b.d12 + a.d12 * b.d22;
		float r20 = a.d20 * b.d00 + a.d21 * b.d10 + a.d22 * b.d20;
		float r21 = a.d20 * b.d01 + a.d21 * b.d11 + a.d22 * b.d21;
		float r22 = a.d20 * b.d02 + a.d21 * b.d12 + a.d22 * b.d22;

		r.d00 = r00;
		r.d01 = r01;
		r.d02 = r02;
		r.d10 = r10;
		r.d11 = r11;
		r.d12 = r12;
		r.d20 = r20;
		r.d21 = r21;
		r.d22 = r22;
	}

	static final private double Quarter = Math.PI / 4.0;
	static final private double Half = Math.PI / 2.0;

	static private double NearSin(double r) {
		double ar = Math.abs(r);
		if (Math.abs(ar - Quarter) >= 1.0e-4) {
			if (Math.abs(ar - Half) >= 1.0e-4) {
				return Math.sin(r);
			} else {
				return 0.0;
			}
		} else {
			return r > 0.0 ? 1.0 : -1.0;
		}
	}

	static private double NearCos(double r) {
		double ar = Math.abs(r);		
		if (Math.abs(ar - Quarter) >= 1.0e-4) {
			if (Math.abs(ar - Half) >= 1.0e-4) {
				return Math.cos(r);
			} else {
				return -1.0;
			}
		} else {
			return 0.0;
		}
	}
}
