// Point3D.h: interface for the CPoint3D class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <math.h>
#define _USE_MATH_DEFINES

struct CPoint3D  
{
public:
	double x, y, z;
	CPoint3D();
	CPoint3D(double x, double y, double z);	
	inline CPoint3D& operator +=(const CPoint3D& pt);
	inline CPoint3D& operator -=(const CPoint3D& pt);
	inline CPoint3D& operator *=(double times);
	inline CPoint3D& operator /=(double times);
	inline CPoint3D operator /(double times) const;
	inline double Len() const;
	inline void Normalize();
	bool operator==(const CPoint3D& other) const
	{
		return x == other.x && y == other.y && z == other.z;
	}
	bool operator<(const CPoint3D& other) const
	{
		if (x < other.x)
			return true;
		else if (x > other.x)
			return false;
		else if (y < other.y)
			return true;
		else if (y > other.y)
			return false;
		else if (z < other.z)
			return true;
		else if (z > other.z)
			return false;
		return false;
	}
};

CPoint3D& CPoint3D::operator +=(const CPoint3D& pt)
{
	x += pt.x;
	y += pt.y;
	z += pt.z;
	return *this;
}

CPoint3D& CPoint3D::operator -=(const CPoint3D& pt)
{
	x -= pt.x;
	y -= pt.y;
	z -= pt.z;
	return *this;
}

CPoint3D& CPoint3D::operator *=(double times)
{
	x *= times;
	y *= times;
	z *= times;
	return *this;
}

CPoint3D& CPoint3D::operator /=(double times)
{
	x /= times;
	y /= times;
	z /= times;
	return *this;
}

CPoint3D CPoint3D::operator /(double times) const
{
	return CPoint3D(x / times, y / times, z / times);
}

double CPoint3D::Len() const
{
	return sqrt(x * x + y * y + z * z);
}

void CPoint3D::Normalize()
{
	double len = Len();
	x /= len;
	y /= len;
	z /= len;
}

CPoint3D operator +(const CPoint3D& pt1, const CPoint3D& pt2);
CPoint3D operator -(const CPoint3D& pt1, const CPoint3D& pt2);
CPoint3D operator *(const CPoint3D& pt, double times);
CPoint3D operator *(double times, const CPoint3D& pt);
CPoint3D operator *(const CPoint3D& pt1, const CPoint3D& pt2);
CPoint3D VectorCross(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
double operator ^(const CPoint3D& pt1, const CPoint3D& pt2);
double GetTriangleArea(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
double AngleBetween(const CPoint3D& pt1, const CPoint3D& pt2);
double AngleBetween(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
void VectorCross(const float* u, const float* v, float * n);
float VectorDot(const float* u, const float* v);
float AngleBetween(const float* u, const float* v);
