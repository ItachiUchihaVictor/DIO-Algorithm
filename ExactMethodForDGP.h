// ExactMethodForDGP.h: interface for the CExactMethodForDGP class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "RichModel.h"
#include <list>
#include <float.h>
using namespace std;
const double LENGTH_EPSILON_CONTROL = 1e-6;


class CExactMethodForDGP  
{
public:
	struct InfoAtVertex
	{
		bool fParentIsPseudoSource;
		char birthTime;		
		int indexOfParent;
		int indexOfRootVertOfParent;
		int level;
		double disUptodate;
		double entryProp;

		InfoAtVertex()
		{
			birthTime = -1;
			disUptodate = DBL_MAX;
		}		
	};
	struct QuoteInfoAtVertex
	{
		char birthTime;
		int indexOfVert;
		double disUptodate;
		double disAstar;
		bool operator<(const QuoteInfoAtVertex& another) const
		{
			return disAstar < another.disAstar;
		}
		QuoteInfoAtVertex(){}
		QuoteInfoAtVertex(char birthTime, int indexOfVert, double disUptodate)
		{
			this->birthTime = birthTime;
			this->indexOfVert = indexOfVert;
			this->disUptodate = disUptodate;
		}
	};
	vector<InfoAtVertex> m_InfoAtVertices;
    int visited_vertices = 0;
    int GetDest(){
        return dest;
    }
    void SetDest(int Y){
        dest = Y;
    }
protected:	
	bool fComputationCompleted;	
	bool fLocked;
	double totalLen;
	int nTotalCurves;

	set<int> indexOfSourceVerts;
	int nCountOfWindows;
	double clockTicks;
	int nMaxLenOfWindowQueue;
	double nMaxLenOfPseudoSources;
	int depthOfResultingTree;
	double NPE;
	double memory;
	double farestDis;
	vector<list<CPoint3D> > m_tableOfResultingPaths;

	string nameOfAlgorithm;
protected:
	void BackTrace(int indexOfVert);
	void BackTraceWithoutStoring(int indexOfVert) const;
private:
    int dest;
public:
	CExactMethodForDGP(const CRichModel& inputModel, const set<int> &indexOfSourceVerts);
	virtual ~CExactMethodForDGP();
	double GetDistanceAt(int v) const {return m_InfoAtVertices[v].disUptodate;}
	inline int GetRootSourceOfVert(int index) const;
	int FindSourceVertex(int indexOfVert, vector<EdgePoint>& resultingPath) const; 
	void PickShortestPaths(int num);
	// virtual void Render() const;
	// virtual void DrawIsolines(int num);
	// void SaveIsolines(int num, const char* filename) const;
	// void SaveIsolines(int num, int order, const char* filename) const;
	// virtual void DrawVoronoi() const;
	virtual void Execute();
	virtual void InitContainers() = 0;
	virtual void BuildSequenceTree() = 0;
	virtual void ClearContainers() = 0;
	virtual void FillExperimentalResults() = 0;
	inline double GetRunTime() const;
	inline double GetMemoryCost() const;
	inline int GetWindowNum() const;
	inline int GetMaxLenOfQue() const;
	inline double GetNPE() const;
	inline int GetDepthOfSequenceTree() const;
	inline string GetAlgorithmName() const;
	inline bool HasBeenCompleted() const;

	const CRichModel& model;
};

double CExactMethodForDGP::GetRunTime() const
{
	return clockTicks;
}

double CExactMethodForDGP::GetMemoryCost() const
{
	return memory;
}

int CExactMethodForDGP::GetWindowNum() const
{
	return nCountOfWindows;
}

int CExactMethodForDGP::GetMaxLenOfQue() const
{
	return nMaxLenOfWindowQueue;
}

int CExactMethodForDGP::GetDepthOfSequenceTree() const
{
	return depthOfResultingTree;
}

double CExactMethodForDGP::GetNPE() const
{
	return NPE;
}

string CExactMethodForDGP::GetAlgorithmName() const
{
	return nameOfAlgorithm;
}

bool CExactMethodForDGP::HasBeenCompleted() const
{
	return fComputationCompleted;
}

int CExactMethodForDGP::GetRootSourceOfVert(int index) const
{
	if (m_InfoAtVertices[index].disUptodate > FLT_MAX)
		return index;

	while (m_InfoAtVertices[index].disUptodate > FLT_EPSILON)
	{
		int indexOfParent = m_InfoAtVertices[index].indexOfParent;
		if (m_InfoAtVertices[index].fParentIsPseudoSource)
		{
			index = indexOfParent;
		}
		else
		{
			index = m_InfoAtVertices[index].indexOfRootVertOfParent;
		}
	}
	return index;
}

