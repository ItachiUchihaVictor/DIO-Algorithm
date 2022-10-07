// BaseModel.h: interface for the CBaseModel class.
//
//////////////////////////////////////////////////////////////////////
#pragma once

#include "Point3D.h"
#include <string>
#include <vector>
#include <map>
#include <set>

#ifndef BASEMODEL_H
#define	BASEMODEL_H

using namespace std;


  
class CBaseModel  
{
public:	
	CBaseModel(const string& filename);
 public:
	struct CFace
	{
		int verts[3];
		CFace(){}
		CFace(int x, int y, int z)
		{
			verts[0] = x;
			verts[1] = y;
			verts[2] = z;
		}
		int& operator[](int index)
		{
			return verts[index];
		}
		int operator[](int index) const
		{
			return verts[index];
		} 
	};
	
protected:	
	void ReadMFile(const string& filename);
	void ReadFile(const string& filename);
	void ReadObjFile(const string& filename);
	void ReadOffFile(const string& filename);
public:
	void AdjustScaleAndComputeNormalsToVerts();
	void SaveMFile(const string& filename) const;
	void SaveOffFile(const string& filename) const;
	void SaveObjFile(const string& filename) const;
	// virtual void Render() const;
	inline int GetNumOfVerts() const;
	inline int GetNumOfFaces() const;
	void LoadModel();
	string GetFileName() const; 
	inline const CPoint3D& Vert(int vertIndex) const;
	inline const CPoint3D& Normal(int vertIndex) const;
	inline const CFace& Face(int faceIndex) const;
	inline bool HasBeenLoad() const;
public:
	vector<CPoint3D> m_Verts;
	vector<CFace> m_Faces;
	vector<CPoint3D> m_NormalsToVerts;
	set<int> m_UselessFaces;
	bool m_fBeLoaded;
	string m_filename;
	CPoint3D m_center;
	double m_scale;
	CPoint3D m_ptUp;
	CPoint3D m_ptDown;
};

int CBaseModel::GetNumOfVerts() const
{
	return (int)m_Verts.size();
}

int CBaseModel::GetNumOfFaces() const
{
	return (int)m_Faces.size();
}

const CPoint3D& CBaseModel::Vert(int vertIndex) const
{
	return m_Verts[vertIndex];
}

const CPoint3D& CBaseModel::Normal(int vertIndex) const
{
	return m_NormalsToVerts[vertIndex];
}

const CBaseModel::CFace& CBaseModel::Face(int faceIndex) const
{
	return m_Faces[faceIndex];
}

bool CBaseModel::HasBeenLoad() const
{
	return m_fBeLoaded;
}

#endif	/* BASEMODEL_H */