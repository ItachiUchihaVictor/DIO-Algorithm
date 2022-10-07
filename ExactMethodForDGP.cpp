// ExactMethodForDGP.cpp: implementation of the CExactMethodForDGP class.
//
//////////////////////////////////////////////////////////////////////
#include "ExactMethodForDGP.h"
//#include <windows.h>
//#include <gl/GL.h>
//#include <gl/GLU.h>
//#pragma comment(lib, "opengl32.lib")
//#pragma comment(lib, "glu32.lib")


// #include <GL/glut.h>

#include <fstream>
#include <time.h>
using namespace std;
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CExactMethodForDGP::CExactMethodForDGP(const CRichModel& inputModel, const set<int> &indexOfSourceVerts) : model(inputModel)
{
	this->indexOfSourceVerts = indexOfSourceVerts;
	nCountOfWindows = 0;
	nMaxLenOfPseudoSources = 0;
	nMaxLenOfWindowQueue = 0;
	depthOfResultingTree = 0;
	totalLen = 0;
	fComputationCompleted = false;
	fLocked = false;
	NPE = 0;
	memory = 0;
	nTotalCurves = 0;
	nameOfAlgorithm = "";
	m_InfoAtVertices.resize(model.GetNumOfVerts());
	memory += double(model.GetNumOfVerts()) * sizeof(InfoAtVertex) / 1024 / 1024;
}

CExactMethodForDGP::~CExactMethodForDGP()
{
}

void CExactMethodForDGP::PickShortestPaths(int num)
{
	if (num >= model.GetNumOfVerts())
		num = model.GetNumOfVerts();
	nTotalCurves = num;
	m_tableOfResultingPaths.clear();
	if (num == 0)
		return;
	if (model.GetNumOfFaces() * num < 4e6)
	{
		if (num >= model.GetNumOfVerts())
		{
			m_tableOfResultingPaths.reserve(model.GetNumOfVerts());
			for (int i = 0; i < model.GetNumOfVerts(); ++i)
			{
				BackTrace(i);
			}
		}
		else
		{
			float step = model.GetNumOfVerts() / float(num);
			//step = max<float>(1, step);
			step = step > 1 ? step : 1;
			m_tableOfResultingPaths.reserve(int(model.GetNumOfVerts() / step) + 1);
			for (float i = FLT_EPSILON; i < model.GetNumOfVerts(); i += step)
			{
				BackTrace(int(i));
			}
		}
	}	
}

void CExactMethodForDGP::BackTrace(int indexOfVert)
{
	if (m_InfoAtVertices[indexOfVert].birthTime == -1)
	{
		assert(model.GetNumOfComponents() != 1 || model.Neigh(indexOfVert).empty());
		return;
	}
	m_tableOfResultingPaths.push_back(list<CPoint3D>());
	vector<int> vertexNodes;
	int index = indexOfVert;
	vertexNodes.push_back(index);
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
		vertexNodes.push_back(index);
	};
	int indexOfSourceVert = index;
	int posOfTable = (int)m_tableOfResultingPaths.size() - 1;
	for (int i = 0; i < (int)vertexNodes.size() - 1; ++i)
	{
		int lastVert = vertexNodes[i];
		CPoint3D pt = model.ComputeShiftPoint(lastVert);
		m_tableOfResultingPaths[posOfTable].push_back(pt);
		
		if (m_InfoAtVertices[lastVert].fParentIsPseudoSource)
		{
			continue;
		}
		int parentEdgeIndex = m_InfoAtVertices[lastVert].indexOfParent;
		int edgeIndex = model.Edge(parentEdgeIndex).indexOfReverseEdge;
		pair<double, double> coord(model.GetNew2DCoordinatesByReversingCurrentEdge(parentEdgeIndex, model.Edge(parentEdgeIndex).coordOfOppositeVert));
		
		double proportion = 1 - m_InfoAtVertices[lastVert].entryProp;
		while (1) 
		{
			CPoint3D pt1 = model.ComputeShiftPoint(model.Edge(edgeIndex).indexOfLeftVert);
			CPoint3D pt2 = model.ComputeShiftPoint(model.Edge(edgeIndex).indexOfRightVert);
			CPoint3D ptIntersection = CRichModel::CombineTwoNormalsTo(pt1, 1 - proportion, pt2, proportion);			
			m_tableOfResultingPaths[posOfTable].push_back(ptIntersection);

			if (model.Edge(edgeIndex).indexOfOppositeVert == vertexNodes[i + 1])
				break;
			double oldProprotion = proportion;
			proportion = model.ProportionOnLeftEdgeByImage(edgeIndex, coord, oldProprotion);
			if (proportion >= -LENGTH_EPSILON_CONTROL && proportion <= 1)
			{
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfLeftEdge;
				//rightLen = disToAngle;				
			}
			else
			{
				proportion = model.ProportionOnRightEdgeByImage(edgeIndex, coord, oldProprotion);
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				//proportion = min<double>(proportion, 1);
				proportion = proportion < 1 ? proportion : 1;
				coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfRightEdge;
			}
		};
	}
	m_tableOfResultingPaths[posOfTable].push_back(model.ComputeShiftPoint(indexOfSourceVert));	
}

void CExactMethodForDGP::BackTraceWithoutStoring(int indexOfVert) const
{
	if (m_InfoAtVertices[indexOfVert].birthTime == -1)
	{
		assert(model.GetNumOfComponents() != 1 || model.Neigh(indexOfVert).empty());
		return;
	}
//	glBegin(GL_LINE_STRIP);
	vector<int> vertexNodes;
	int index = indexOfVert;
	vertexNodes.push_back(index);
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
		vertexNodes.push_back(index);
	};
	int indexOfSourceVert = index;

	for (int i = 0; i < (int)vertexNodes.size() - 1; ++i)
	{
		int lastVert = vertexNodes[i];
		CPoint3D pt = model.ComputeShiftPoint(lastVert);
//		glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
		if (m_InfoAtVertices[lastVert].fParentIsPseudoSource)
		{
			continue;
		}
		int parentEdgeIndex = m_InfoAtVertices[lastVert].indexOfParent;
		int edgeIndex = model.Edge(parentEdgeIndex).indexOfReverseEdge;
		pair<double, double> coord(model.GetNew2DCoordinatesByReversingCurrentEdge(parentEdgeIndex, model.Edge(parentEdgeIndex).coordOfOppositeVert));
		
		double proportion = 1 - m_InfoAtVertices[lastVert].entryProp;
		while (1) 
		{
			CPoint3D pt1 = model.ComputeShiftPoint(model.Edge(edgeIndex).indexOfLeftVert);
			CPoint3D pt2 = model.ComputeShiftPoint(model.Edge(edgeIndex).indexOfRightVert);
			CPoint3D ptIntersection = CRichModel::CombineTwoNormalsTo(pt1, 1 - proportion, pt2, proportion);
//			glVertex3f((float)ptIntersection.x, (float)ptIntersection.y, (float)ptIntersection.z);
			if (model.Edge(edgeIndex).indexOfOppositeVert == vertexNodes[i + 1])
				break;
			double oldProprotion = proportion;
			proportion = model.ProportionOnLeftEdgeByImage(edgeIndex, coord, oldProprotion);
			if (proportion >= -LENGTH_EPSILON_CONTROL && proportion <= 1)
			{
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfLeftEdge;
				//rightLen = disToAngle;				
			}
			else
			{
				proportion = model.ProportionOnRightEdgeByImage(edgeIndex, coord, oldProprotion);
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				//proportion = min<double>(proportion, 1);
				proportion = proportion < 1 ? proportion : 1;
				coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfRightEdge;
			}
		};
	}
	CPoint3D pt = model.ComputeShiftPoint(indexOfSourceVert);	
//	glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
//	glEnd();
}

int CExactMethodForDGP::FindSourceVertex(int indexOfVert, vector<EdgePoint>& resultingPath) const
{
	resultingPath.clear();

	if (m_InfoAtVertices[indexOfVert].birthTime == -1 || m_InfoAtVertices[indexOfVert].disUptodate > FLT_MAX)
	{
		assert(model.GetNumOfComponents() != 1 || model.Neigh(indexOfVert).empty());
		return -1;
	}
	vector<int> vertexNodes;
	int index = indexOfVert;
	vertexNodes.push_back(index);
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
		vertexNodes.push_back(index);
	};
	int indexOfSourceVert = index;

	for (int i = 0; i < (int)vertexNodes.size() - 1; ++i)
	{
		int lastVert = vertexNodes[i];
		//if (lastVert != indexOfVert)
		resultingPath.push_back(EdgePoint(lastVert));
		if (m_InfoAtVertices[lastVert].fParentIsPseudoSource)
		{
			continue;
		}
		int parentEdgeIndex = m_InfoAtVertices[lastVert].indexOfParent;
		int edgeIndex = model.Edge(parentEdgeIndex).indexOfReverseEdge;
		pair<double, double> coord(model.GetNew2DCoordinatesByReversingCurrentEdge(parentEdgeIndex, model.Edge(parentEdgeIndex).coordOfOppositeVert));

		double proportion = 1 - m_InfoAtVertices[lastVert].entryProp;
		while (1) 
		{
			resultingPath.push_back(EdgePoint(edgeIndex, proportion));
			if (model.Edge(edgeIndex).indexOfOppositeVert == vertexNodes[i + 1])
				break;
			double oldProprotion = proportion;
			proportion = model.ProportionOnLeftEdgeByImage(edgeIndex, coord, oldProprotion);
			if (model.Edge(edgeIndex).indexOfLeftEdge == -1 || model.Edge(edgeIndex).indexOfRightEdge == -1)
			{
				break;
			}

			if (proportion >= -LENGTH_EPSILON_CONTROL && proportion <= 1)
			{
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfLeftEdge;
				//rightLen = disToAngle;			
			}
			else
			{
				proportion = model.ProportionOnRightEdgeByImage(edgeIndex, coord, oldProprotion);
				//proportion = max<double>(proportion, 0);
				proportion = proportion > 0 ? proportion : 0;
				//proportion = min<double>(proportion, 1);
				proportion = proportion < 1 ? proportion : 1;
				coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(edgeIndex, coord);
				edgeIndex = model.Edge(edgeIndex).indexOfRightEdge;
			}
		};
	}
	resultingPath.push_back(EdgePoint(indexOfSourceVert));
	return indexOfSourceVert;
}

/*
void CExactMethodForDGP::Render() const
{
//	glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	srand(0);
	if (!m_tableOfResultingPaths.empty())
	{
		for (int i = 0; i < (int)m_tableOfResultingPaths.size(); ++i)
		{
			double r = rand() / double(RAND_MAX);
			double g = rand() / double(RAND_MAX);
			double b = rand() / double(RAND_MAX);
			glColor3f((float)r, (float)g, (float)b);

			glBegin(GL_LINE_STRIP);
			for (list<CPoint3D>::const_iterator it = m_tableOfResultingPaths[i].begin();
				it != m_tableOfResultingPaths[i].end(); ++it)
			{
				glVertex3f((float)it->x, (float)it->y, (float)it->z);
			}
			glEnd();		
		}
	}
	else if (nTotalCurves != 0)
	{
		float step = model.GetNumOfVerts() / float(nTotalCurves);
		step = max<float>(1, step);
		for (float i = FLT_EPSILON; i < model.GetNumOfVerts(); i += step)
		{
			double r = rand() / double(RAND_MAX);
			double g = rand() / double(RAND_MAX);
			double b = rand() / double(RAND_MAX);
			glColor3f((float)r, (float)g, (float)b);
			BackTraceWithoutStoring(int(i));
		}
	}
	glPopAttrib();	
}
 */

void CExactMethodForDGP::Execute()
{
	if (fComputationCompleted)
		return;
	if (!fLocked)
	{
		fLocked = true;
		nCountOfWindows = 0;	
		nMaxLenOfWindowQueue = 0;
		depthOfResultingTree = 0;
		InitContainers();
		clockTicks = clock();	 // 0; //GetTickCount();	
		BuildSequenceTree();
	  clockTicks = clock() - clockTicks;
		FillExperimentalResults();
		ClearContainers();
		
		fComputationCompleted = true;
		fLocked = false;
	}
}

