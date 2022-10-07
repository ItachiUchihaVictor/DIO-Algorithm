// RichModel.h: interface for the CRichModel class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "BaseModel.h"
const double RateOfNormalShift = 5e-3;
const double ToleranceOfConvexAngle = 5e-3;
#include <cassert>
using namespace std;

class CRichModel : virtual public CBaseModel 
{
public:
	struct CEdge
	{
		int indexOfLeftVert;
		int indexOfRightVert;
		int indexOfOppositeVert;
		int indexOfLeftEdge;
		int indexOfRightEdge;
		int indexOfReverseEdge;
		int indexOfFrontFace;
		double length;
		pair<double, double> coordOfOppositeVert;
		// |unitX   -unitY|
		// |unitY    unitX|
		pair<double, double> matrixRotatedToLeftEdge;
		pair<double, double> matrixRotatedToRightEdge;
		CEdge()
		{		
			indexOfOppositeVert = -1;	//key	
			indexOfLeftEdge = -1;
			indexOfRightEdge = -1;
			indexOfFrontFace = -1;
		}
	};

protected:
	void CreateEdgesFromVertsAndFaces();
	void CollectAndArrangeNeighs();
	void ComputeAnglesAroundVerts();
	void ComputePlanarCoordsOfIncidentVertForEdges();
	void ComputeNumOfHoles();
	void ComputeNumOfComponents();
public:
	CRichModel();
	CRichModel(const string &filename);
	void Preprocess();	
	void ClearEdges() {
		fLocked = false; fBePreprocessed = false; m_nBoundries = 0;
		m_nIsolatedVerts = 0; 
		m_NeighsAndAngles.clear(); m_FlagsForCheckingConvexVerts.clear(); m_Edges.clear();
	}
	inline int GetSubindexToVert(int root, int neigh) const;
	inline const CEdge& Edge(int edgeIndex) const;	
	inline const vector<pair<int, double> >& Neigh(int root) const;	
	inline double AngleSum(int vertIndex) const;
	//inline double Curvature(int vertIndex) const;
	inline double ProportionOnEdgeByImage(int edgeIndex, const pair<double, double> &coord) const;
	inline double ProportionOnLeftEdgeByImage(int edgeIndex, const pair<double, double> &coord, double proportion) const;
	inline double ProportionOnRightEdgeByImage(int edgeIndex, const pair<double, double> &coord, double proportion) const;
	inline double ProportionOnEdgeByImage(int edgeIndex, double x1, double y1, double x2, double y2) const;
	inline pair<double, double> GetNew2DCoordinatesByRotatingAroundLeftChildEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const;	
	inline pair<double, double> GetNew2DCoordinatesByRotatingAroundRightChildEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const;
	inline pair<double, double> GetNew2DCoordinatesByReversingCurrentEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const;
	inline double DistanceToIncidentAngle(int edgeIndex, const pair<double, double>& coord) const;
	inline int GetNumOfEdges() const;
	inline int GetNumOfValidDirectedEdges() const;
	inline int GetNumOfTotalUndirectedEdges() const;
	inline int GetNumOfGenera() const;
	inline int GetNumOfIsolated() const;
	inline int GetNumOfComponents() const;
	inline int GetNumOfBoundries() const;
	inline bool IsConvexVert(int index) const;
	inline bool isBoundaryVert(int index) const;
	inline bool IsClosedModel() const;
	inline bool IsExtremeEdge(int edgeIndex) const;
	inline bool IsStartEdge(int edgeIndex) const;	
	inline bool HasBeenProcessed() const;
	//inline int GetFirstEdgeIndex(int faceIndex) const;
	//inline int GetSecondEdgeIndex(int faceIndex) const;
	//inline int GetThirdEdgeIndex(int faceIndex) const;
	//inline int GetEdgeIndexFromFace(int faceIndex, int subIndex) const;
	inline int GetEdgeIndexFromTwoVertices(int leftVert, int rightVert) const;
	inline CPoint3D ComputeShiftPoint(int indexOfVert) const;
	inline CPoint3D ComputeShiftPoint(int indexOfVert, double epsilon) const;
	static CPoint3D CombinePointAndNormalTo(const CPoint3D& pt, const CPoint3D& normal);
	static CPoint3D CombineTwoNormalsTo(const CPoint3D& pt1, double coef1, const CPoint3D& pt2, double coef2);	
protected:
	bool fLocked;
	bool fBePreprocessed;
	int m_nBoundries;
	int m_nIsolatedVerts;
	int m_nComponents;
public:
	vector<CEdge> m_Edges;
	set<int> m_UselessEdges;
	vector<vector<pair<int, double> > > m_NeighsAndAngles;	
	vector<bool> m_FlagsForCheckingConvexVerts;
};

int CRichModel::GetNumOfValidDirectedEdges() const
{
	return (int)m_Faces.size() * 3;
}

int CRichModel::GetNumOfTotalUndirectedEdges() const
{
	return (int)m_Edges.size() / 2;
}

int CRichModel::GetNumOfGenera() const
{
	return int(GetNumOfTotalUndirectedEdges() - (GetNumOfVerts() - m_nIsolatedVerts) - GetNumOfFaces() - GetNumOfBoundries()) / 2 + 1;
}

int CRichModel::GetNumOfComponents() const
{
	return m_nComponents;
}

int CRichModel::GetNumOfBoundries() const
{
	return m_nBoundries;
}

bool CRichModel::IsClosedModel() const
{
	return GetNumOfValidDirectedEdges() ==  GetNumOfEdges();
}

int CRichModel::GetNumOfIsolated() const
{
	return m_nIsolatedVerts;
}

int CRichModel::GetNumOfEdges() const
{
	return (int)m_Edges.size();
}

bool CRichModel::isBoundaryVert(int index) const
{
	return IsStartEdge(Neigh(index).front().first);
}

bool CRichModel::IsConvexVert(int index) const
{
	return m_FlagsForCheckingConvexVerts[index];
}

bool CRichModel::IsExtremeEdge(int edgeIndex) const
{
	return Edge(edgeIndex).indexOfOppositeVert == -1;
}

bool CRichModel::IsStartEdge(int edgeIndex) const
{
	return Edge(Edge(edgeIndex).indexOfReverseEdge).indexOfOppositeVert == -1;
}

const CRichModel::CEdge& CRichModel::Edge(int edgeIndex) const
{
	return m_Edges[edgeIndex];
}

const vector<pair<int, double> >& CRichModel::Neigh(int root) const
{
	return m_NeighsAndAngles[root];
}

double CRichModel::ProportionOnEdgeByImage(int edgeIndex, const pair<double, double>& coord) const
{
	double res = Edge(edgeIndex).coordOfOppositeVert.first * coord.second - Edge(edgeIndex).coordOfOppositeVert.second * coord.first;
	return res / ((coord.second - Edge(edgeIndex).coordOfOppositeVert.second) * Edge(edgeIndex).length);
}

double CRichModel::ProportionOnEdgeByImage(int edgeIndex, double x1, double y1, double x2, double y2) const
{
	double res = x1 * y2 - x2 * y1;
	return res / ((y2 - y1) * Edge(edgeIndex).length);
}

double CRichModel::ProportionOnLeftEdgeByImage(int edgeIndex, const pair<double, double> &coord, double proportion) const
{
	double xBalance = proportion * Edge(edgeIndex).length;
	double res = Edge(edgeIndex).coordOfOppositeVert.first * coord.second - Edge(edgeIndex).coordOfOppositeVert.second * (coord.first - xBalance);
	return xBalance * coord.second / res;
}

double CRichModel::ProportionOnRightEdgeByImage(int edgeIndex, const pair<double, double> &coord, double proportion) const
{
	double part1 = Edge(edgeIndex).length * coord.second;
	double part2 = proportion * Edge(edgeIndex).length * Edge(edgeIndex).coordOfOppositeVert.second;
	double part3 = Edge(edgeIndex).coordOfOppositeVert.second * coord.first - Edge(edgeIndex).coordOfOppositeVert.first * coord.second;	
	return (part3 + proportion * part1 - part2) / (part3 + part1 - part2);
}

pair<double, double> CRichModel::GetNew2DCoordinatesByRotatingAroundLeftChildEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const
{
	return make_pair(Edge(edgeIndex).matrixRotatedToLeftEdge.first * input2DCoordinates.first - Edge(edgeIndex).matrixRotatedToLeftEdge.second * input2DCoordinates.second,
		Edge(edgeIndex).matrixRotatedToLeftEdge.second * input2DCoordinates.first + Edge(edgeIndex).matrixRotatedToLeftEdge.first * input2DCoordinates.second);
}

pair<double, double> CRichModel::GetNew2DCoordinatesByRotatingAroundRightChildEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const
{
	int reverseEdge = Edge(Edge(edgeIndex).indexOfRightEdge).indexOfReverseEdge;
	pair<double, double> coordOfLeftEnd = GetNew2DCoordinatesByReversingCurrentEdge(reverseEdge, Edge(reverseEdge).coordOfOppositeVert);	
	return make_pair(Edge(edgeIndex).matrixRotatedToRightEdge.first * input2DCoordinates.first - Edge(edgeIndex).matrixRotatedToRightEdge.second * input2DCoordinates.second + coordOfLeftEnd.first,
		Edge(edgeIndex).matrixRotatedToRightEdge.second * input2DCoordinates.first + Edge(edgeIndex).matrixRotatedToRightEdge.first * input2DCoordinates.second + coordOfLeftEnd.second);
}

pair<double, double> CRichModel::GetNew2DCoordinatesByReversingCurrentEdge(int edgeIndex, const pair<double, double>& input2DCoordinates) const
{
	return make_pair(Edge(edgeIndex).length - input2DCoordinates.first, - input2DCoordinates.second);
}

bool CRichModel::HasBeenProcessed() const
{
	return fBePreprocessed;
}

int CRichModel::GetSubindexToVert(int root, int neigh) const
{
	for (int i = 0; i < (int)Neigh(root).size(); ++i)
	{
		if (Edge(Neigh(root)[i].first).indexOfRightVert == neigh)
			return i;
	}
	return -1;
}

CPoint3D CRichModel::ComputeShiftPoint(int indexOfVert) const
{
	return Vert(indexOfVert) + Normal(indexOfVert) * RateOfNormalShift / m_scale;
}

CPoint3D CRichModel::ComputeShiftPoint(int indexOfVert, double epsilon) const
{
	return Vert(indexOfVert) +  Normal(indexOfVert) * epsilon;
}

double CRichModel::AngleSum(int vertIndex) const
{
	double angleSum(0);
	for (int j = 0; j < (int)m_NeighsAndAngles[vertIndex].size(); ++j)
	{		
		angleSum += m_NeighsAndAngles[vertIndex][j].second;			
	}
	return angleSum;
}

double CRichModel::DistanceToIncidentAngle(int edgeIndex, const pair<double, double>& coord) const
{
	double detaX = coord.first - Edge(edgeIndex).coordOfOppositeVert.first;
	double detaY = coord.second - Edge(edgeIndex).coordOfOppositeVert.second;
	return sqrt(detaX * detaX + detaY * detaY);
}

//int CRichModel::GetFirstEdgeIndex(int faceIndex) const
//{
//	int root = m_Faces[faceIndex][0];
//	int subIndex = GetSubindexToVert(root, m_Faces[faceIndex][1]);
//	return Neigh(root)[subIndex].first;
//}
//int CRichModel::GetSecondEdgeIndex(int faceIndex) const
//{
//	int root = m_Faces[faceIndex][1];
//	int subIndex = GetSubindexToVert(root, m_Faces[faceIndex][2]);
//	return Neigh(root)[subIndex].first;
//}
//int CRichModel::GetThirdEdgeIndex(int faceIndex) const
//{
//	int root = m_Faces[faceIndex][2];
//	int subIndex = GetSubindexToVert(root, m_Faces[faceIndex][0]);
//	return Neigh(root)[subIndex].first;
//}
//int CRichModel::GetEdgeIndexFromFace(int faceIndex, int subIndex) const
//{
//	if (subIndex == 0)
//	{
//		int edgeIndex = GetFirstEdgeIndex(faceIndex);
//		assert (Edge(edgeIndex).indexOfFrontFace == faceIndex);
//		return edgeIndex;
//	}
//	if (subIndex == 1)
//	{
//		int edgeIndex = GetSecondEdgeIndex(faceIndex);
//		assert (Edge(edgeIndex).indexOfFrontFace == faceIndex);
//		return edgeIndex;
//	}
//	if (subIndex == 2)
//	{
//		int edgeIndex = GetThirdEdgeIndex(faceIndex);
//		assert (Edge(edgeIndex).indexOfFrontFace == faceIndex);
//		return edgeIndex;
//	}
//	assert(false);
//	return -1;
//}
//
int CRichModel::GetEdgeIndexFromTwoVertices(int leftVert, int rightVert) const
{
	int subIndex = GetSubindexToVert(leftVert, rightVert);
	assert (subIndex != -1);
	return Neigh(leftVert)[subIndex].first;
}


struct EdgePoint
{
	bool isVertex;
	int index;
	double proportion; //[0 --> left endpoint; 1 --> right endpoint]
	EdgePoint()
	{
	}
	EdgePoint(int index) : index(index), isVertex(true){}
	EdgePoint(int index, double proportion) : index(index), proportion(proportion), isVertex(false) {}
	EdgePoint(const CRichModel& model, int leftVert, int rightVert, double proportion) : proportion(proportion), isVertex(false) 
	{
		index = model.GetEdgeIndexFromTwoVertices(leftVert, rightVert);
	}
	CPoint3D Get3DPoint(const CRichModel& model)
	{
		if (isVertex)
			return model.Vert(index);
		return (1 - proportion) * model.Vert(model.Edge(index).indexOfLeftVert)
			+ proportion * model.Vert(model.Edge(index).indexOfRightVert);
	}
	bool operator <(const EdgePoint& other) const
	{
		if (isVertex == false && other.isVertex == true)
			return true;
		if (isVertex == true && other.isVertex == false)
			return false;
		if (index < other.index)
			return true;
		if (index > other.index)
			return false;
		if (proportion < other.proportion)
			return true;
		if (proportion > other.proportion)
			return false;
		return false;
	}
	int GetLeftVertIndex(const CRichModel& model) const
	{
		return model.Edge(index).indexOfLeftVert;
	}
	int GetRightVertIndex(const CRichModel& model) const
	{
		return model.Edge(index).indexOfRightVert;
	}
};

