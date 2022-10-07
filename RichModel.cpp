#include "RichModel.h"
#include <queue>
#include <cmath>
#include <iostream>
#define _USE_MATH_DEFINES
using namespace std;

CRichModel::CRichModel(const string& filename) : CBaseModel(filename)
{
	fBePreprocessed = false;
	fLocked = false;
}

CRichModel::CRichModel() : CBaseModel("")
{
	fBePreprocessed = false;
	fLocked = false;
}


void CRichModel::CreateEdgesFromVertsAndFaces()
{
	m_Edges.reserve(2 * (GetNumOfVerts() + GetNumOfFaces() - 2));
	map<pair<int, int>, int> pondOfUndeterminedEdges;
	int szFaces = GetNumOfFaces();
	for (int i = 0; i < szFaces; ++i)
	{		
		int threeIndices[3];
		for (int j = 0; j < 3; ++j)
		{
			int post = (j + 1) % 3;
			int pre = (j + 2) % 3;
			
			int leftVert = Face(i)[pre];
			int rightVert = Face(i)[j];

			map<pair<int, int>, int>::const_iterator it = pondOfUndeterminedEdges.find(make_pair(leftVert, rightVert));
			if (it != pondOfUndeterminedEdges.end())
			{
				int posInEdgeList = it->second;
				if (m_Edges[posInEdgeList].indexOfOppositeVert != -1)
				{
          cout << "WTF\n";
					throw "Repeated edges!";
				}
				threeIndices[j] = posInEdgeList;
				m_Edges[posInEdgeList].indexOfOppositeVert = Face(i)[post];
				m_Edges[posInEdgeList].indexOfFrontFace = i;				
			}
			else
			{
				CEdge edge;
				edge.indexOfLeftVert = leftVert;
				edge.indexOfRightVert = rightVert;
				edge.indexOfFrontFace = i;
				edge.indexOfOppositeVert = Face(i)[post];
				edge.indexOfReverseEdge = (int)m_Edges.size() + 1;
				edge.length = (Vert(leftVert) - Vert(rightVert)).Len();
				m_Edges.push_back(edge);
				pondOfUndeterminedEdges[make_pair(leftVert, rightVert)] = threeIndices[j] = (int)m_Edges.size() - 1;

				edge.indexOfLeftVert = rightVert;
				edge.indexOfRightVert = leftVert;
				edge.indexOfReverseEdge = (int)m_Edges.size() - 1;
				edge.indexOfOppositeVert = -1;
				edge.indexOfFrontFace = -1;
				m_Edges.push_back(edge);
				pondOfUndeterminedEdges[make_pair(rightVert, leftVert)] = (int)m_Edges.size() - 1;
			}
		}
		for (int j = 0; j < 3; ++j)
		{
			m_Edges[threeIndices[j]].indexOfLeftEdge = Edge(threeIndices[(j + 2) % 3]).indexOfReverseEdge;
			m_Edges[threeIndices[j]].indexOfRightEdge = Edge(threeIndices[(j + 1) % 3]).indexOfReverseEdge;
		}
	}
	// m_Edges.swap(vector<CEdge>(m_Edges));
  vector<CEdge>(m_Edges).swap(m_Edges);
}

void CRichModel::CollectAndArrangeNeighs()
{
	m_nIsolatedVerts = 0;
	vector<int> sequenceOfDegrees(GetNumOfVerts(), 0);	
	m_NeighsAndAngles.resize(GetNumOfVerts());
	for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
	{
		m_NeighsAndAngles[i].resize(1, make_pair(-1, 0));
	}
	for (int i = 0; i < (int)GetNumOfEdges(); ++i)
	{
		const CEdge& edge = Edge(i);
		++sequenceOfDegrees[edge.indexOfLeftVert];
		int &indexOfStartEdge = m_NeighsAndAngles[edge.indexOfLeftVert][0].first;
		if (indexOfStartEdge == -1 || !IsStartEdge(indexOfStartEdge))
		{
			indexOfStartEdge = i;
		}
		else if (IsStartEdge(i))
		{
			m_NeighsAndAngles[edge.indexOfLeftVert].push_back(make_pair(i, 0));
		}
	}
	for (int i = 0; i < GetNumOfVerts(); ++i)
	{
		if (m_NeighsAndAngles[i][0].first == -1)
		{
			m_NeighsAndAngles[i].clear();
			m_nIsolatedVerts++;	
			continue;
		}
		vector<int> startEdges;
		for (int j = 0; j < (int)Neigh(i).size(); ++j)
		{
			startEdges.push_back(Neigh(i)[j].first);
		}	
		m_NeighsAndAngles[i].resize(sequenceOfDegrees[i], make_pair(0, 0));
		int num(0);
		for (int j = 0; j < (int)startEdges.size(); ++j)
		{
			int curEdge = startEdges[j];			
			while (1)
			{
				m_NeighsAndAngles[i][num].first = curEdge;
				++num;
				if (num >= sequenceOfDegrees[i])
					break;
				if (IsExtremeEdge(curEdge))
					break;
				curEdge = Edge(curEdge).indexOfLeftEdge;
				if (curEdge == startEdges[j])
				{
					break;
				}
			}
		}
		if (num != sequenceOfDegrees[i])
		{
			throw "Complex vertices";
		}
	}
}

void CRichModel::ComputeAnglesAroundVerts()
{	
	m_FlagsForCheckingConvexVerts.resize(GetNumOfVerts());
	for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
	{
		m_NeighsAndAngles[i].resize(Neigh(i).size());
	}
	for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
	{
		double angleSum(0);
		for (int j = 0; j < (int)m_NeighsAndAngles[i].size(); ++j)
		{
			if (IsExtremeEdge(Neigh(i)[j].first))
				m_NeighsAndAngles[i][j].second = 2 * M_PI + 0.1;
			else
			{
				int next = j + 1;
				if (next >= (int)m_NeighsAndAngles[i].size())
				{
					next = 0;
				}
				double l = Edge(Neigh(i)[j].first).length;
				double r = Edge(Neigh(i)[next].first).length;
				double b = Edge(Edge(Neigh(i)[j].first).indexOfRightEdge).length;				
				m_NeighsAndAngles[i][j].second = acos((l * l + r * r - b * b) / (2 * l * r));
			}
			angleSum += m_NeighsAndAngles[i][j].second;			
		}
		m_FlagsForCheckingConvexVerts[i] = (angleSum < 2 * M_PI - ToleranceOfConvexAngle);
	}
}

void CRichModel::ComputePlanarCoordsOfIncidentVertForEdges()
{
	for (int i = 0; i < GetNumOfEdges(); ++i)
	{
		if (IsExtremeEdge(i))
			continue;
		double bottom = Edge(i).length;
		double leftLen = Edge(Edge(i).indexOfLeftEdge).length;
		double squareOfLeftLen = leftLen * leftLen;
		double rightLen = Edge(Edge(i).indexOfRightEdge).length;
		double x = (squareOfLeftLen - rightLen * rightLen) / bottom + bottom;
		x /= 2.0;
		m_Edges[i].coordOfOppositeVert = make_pair(x, sqrt(max(0.0, squareOfLeftLen - x * x)));
	}
	for (int i = 0; i < GetNumOfEdges(); ++i)
	{
		if (IsExtremeEdge(i))
			continue;
		{
			int reverseEdge = m_Edges[m_Edges[i].indexOfLeftEdge].indexOfReverseEdge;
			pair<double, double> coord = GetNew2DCoordinatesByReversingCurrentEdge(reverseEdge, m_Edges[reverseEdge].coordOfOppositeVert);
			double scale = abs(coord.first) + abs(coord.second);
			coord.first /= scale;
			coord.second /= scale;
			double len = sqrt(coord.first * coord.first + coord.second * coord.second);
			m_Edges[i].matrixRotatedToLeftEdge = make_pair(coord.first / len, coord.second / len);
		}
		{
			int reverseEdge = m_Edges[m_Edges[i].indexOfRightEdge].indexOfReverseEdge;
			double rightX = m_Edges[reverseEdge].length;
			double rightY = 0;
			double leftX = m_Edges[reverseEdge].length - m_Edges[reverseEdge].coordOfOppositeVert.first;
			double leftY = -m_Edges[reverseEdge].coordOfOppositeVert.second;

			double detaX = rightX - leftX;
			double detaY = rightY - leftY;
			double scale = abs(detaX) + abs(detaY);
			detaX /= scale;
			detaY /= scale;
			double len = sqrt(detaX * detaX + detaY * detaY);
			m_Edges[i].matrixRotatedToRightEdge = make_pair(detaX / len, detaY / len);
		}
	}
}
/*
void CRichModel::Preprocess()
{
	ClearEdges();
	if (fBePreprocessed)
		return;	
	if (!m_fBeLoaded)
	{
		LoadModel();
	}

	if (!fLocked)
	{
		fLocked = true;
		CreateEdgesFromVertsAndFaces();
		CollectAndArrangeNeighs();	
		ComputeNumOfHoles();
		ComputeAnglesAroundVerts();
		ComputePlanarCoordsOfIncidentVertForEdges();
		ComputeNumOfComponents();
		fBePreprocessed = true;
		fLocked = false; 		
	}
}
*/
void CRichModel::Preprocess()
{
  //std::cout << "[in Preprocess:]" << std::endl;
	if (fBePreprocessed)
		return;	
	if (!m_fBeLoaded)
	{
		LoadModel();
	}
  //std::cout << "[in Preprocess:] After LoadModel()" << endl;
  

	if (!fLocked)
	{
    //std::cout << "[in Preprocess:] flocked is false" << endl;

		fLocked = true;
    
    //std::cout << "[in Preprocess:] Before CreateEdgesFromVertsAndFaces()" << endl;
		CreateEdgesFromVertsAndFaces();
    //std::cout << "[in Preprocess:] After CreateEdgesFromVertsAndFaces()" << endl;

		CollectAndArrangeNeighs();	
    //std::cout << "[in Preprocess:] After CollectAndArrangeNeighs()" << endl;
    
		ComputeNumOfHoles();
    //std::cout << "[in Preprocess:] After ComputeNumOfHoles()" << endl;
    
		ComputeAnglesAroundVerts();
    //std::cout << "[in Preprocess:] After ComputeAnglesAroundVerts()" << endl;
    
		ComputePlanarCoordsOfIncidentVertForEdges();
    //std::cout << "[in Preprocess:] After omputePlanarCoordsOfIncidentVertForEdges()" << endl;
    
		fBePreprocessed = true;
		fLocked = false; 
	}
}

CPoint3D CRichModel::CombinePointAndNormalTo(const CPoint3D& pt, const CPoint3D& normal)
{
	return pt + normal * RateOfNormalShift;
}

CPoint3D CRichModel::CombineTwoNormalsTo(const CPoint3D& pt1, double coef1, const CPoint3D& pt2, double coef2)
{
	return coef1 * pt1 + coef2 * pt2;
}

void CRichModel::ComputeNumOfHoles()
{
	m_nBoundries = 0;
	if (IsClosedModel())
	{
		return;
	}
	set<int> extremeEdges;
	for (int i = 0; i < (int)m_Edges.size(); ++i)
	{
		if (m_Edges[i].indexOfOppositeVert != -1)
			continue;
		extremeEdges.insert(i);
	}		

	while (!extremeEdges.empty())
	{
		++m_nBoundries;
		int firstEdge = *extremeEdges.begin();
		int edge = firstEdge;
		do
		{			
			extremeEdges.erase(edge);
			int root = Edge(edge).indexOfRightVert;
			int index = GetSubindexToVert(root, Edge(edge).indexOfLeftVert);
			edge  = Neigh(root)[(index - 1 + (int)Neigh(root).size()) % (int)Neigh(root).size()].first;		
		} while (edge != firstEdge && !extremeEdges.empty());
	}
}

void CRichModel::ComputeNumOfComponents()
{
	m_nComponents = 0;
	vector<bool> flags(GetNumOfVerts(), false);
	int cnt(0);
	while (cnt < GetNumOfVerts())
	{
		int v;
		for (int i = 0; i < (int)flags.size(); ++i)
		{
			if (!flags[i]) 
			{
				v = i;
				break;
			}
		}
		queue<int> Que;
		Que.push(v);
		while (!Que.empty())
		{
			int v = Que.front();
			Que.pop();
			if (flags[v])
				continue;
			flags[v] = true;
			cnt++;
			for (int i = 0; i < (int)Neigh(v).size(); ++i)
			{
				if (!flags[Edge(Neigh(v)[i].first).indexOfRightVert])
				{
					Que.push(Edge(Neigh(v)[i].first).indexOfRightVert);
				}
			}
		}
		m_nComponents++;
	}
}
