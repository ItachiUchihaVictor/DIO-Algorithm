// ImprovedCHWithEdgeValve.cpp: implementation of the CImprovedCHWithEdgeValve class.
//

#include "ImprovedCHWithEdgeValve.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


CImprovedCHWithEdgeValve::CImprovedCHWithEdgeValve(const CRichModel& inputModel, const set<int> &indexOfSourceVerts) : CPreviousCH(inputModel, indexOfSourceVerts)
{
	nameOfAlgorithm = "ICH1";
}
void CImprovedCHWithEdgeValve::BuildSequenceTree()
{
	ComputeChildrenOfSource();
	bool fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();
	while (depthOfResultingTree < model.GetNumOfFaces() && !(m_QueueForPseudoSources.empty() && m_QueueForWindows.empty()))
	{		
		if ((int)m_QueueForWindows.size() > nMaxLenOfWindowQueue)
			nMaxLenOfWindowQueue = (int)m_QueueForWindows.size();
		if (m_QueueForPseudoSources.size() > nMaxLenOfPseudoSources)
			nMaxLenOfPseudoSources = (int)m_QueueForPseudoSources.size();
		if (fFromQueueOfPseudoSources) //pseudosource
		{				
			int indexOfVert = (*m_QueueForPseudoSources.begin()).indexOfVert;
			m_QueueForPseudoSources.erase(m_QueueForPseudoSources.begin());//pop();			
			ComputeChildrenOfPseudoSource(indexOfVert);			
		}
		else			
		{
			QuoteWindow quoteW = *m_QueueForWindows.begin();
			m_QueueForWindows.erase(m_QueueForWindows.begin());//pop();
			ComputeChildrenOfWindow(quoteW);		
			delete quoteW.pWindow;
		}
		fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();
        //if(GetDistanceAt(GetDest())<DBL_MAX) break;
        if(m_InfoAtVertices[GetDest()].birthTime>-1) break;
	}
}

CImprovedCHWithEdgeValve::~CImprovedCHWithEdgeValve()
{

}

bool CImprovedCHWithEdgeValve::CheckValidityOfWindow(Window& w)
{
	if (w.fDirectParenIsPseudoSource)
		return true;
	const CRichModel::CEdge& edge = model.Edge(w.indexOfCurEdge);
	int leftVert = edge.indexOfLeftVert;
	double detaX = w.coordOfPseudoSource.first - w.proportions[1] * edge.length;
	double rightLen = sqrt(detaX * detaX + w.coordOfPseudoSource.second * w.coordOfPseudoSource.second);
	if (m_InfoAtVertices[leftVert].disUptodate < 10000  / model.m_scale && m_InfoAtVertices[leftVert].disUptodate + w.proportions[1] * edge.length
		< w.disToRoot + rightLen)
	{
		return false;
	}
	int rightVert = edge.indexOfRightVert;
	detaX = w.coordOfPseudoSource.first - w.proportions[0] * edge.length;
	double leftLen = sqrt(detaX * detaX + w.coordOfPseudoSource.second * w.coordOfPseudoSource.second);
	if (m_InfoAtVertices[rightVert].disUptodate < 10000  / model.m_scale && m_InfoAtVertices[rightVert].disUptodate + (1 - w.proportions[0]) * edge.length
		< w.disToRoot + leftLen)
	{
		return false;
	}
	const CRichModel::CEdge& oppositeEdge = model.Edge(edge.indexOfReverseEdge);
	double xOfVert = edge.length - oppositeEdge.coordOfOppositeVert.first;
	double yOfVert = -oppositeEdge.coordOfOppositeVert.second;
	if (m_InfoAtVertices[oppositeEdge.indexOfOppositeVert].disUptodate < 10000  / model.m_scale)	
	{
		if (w.fDirectParentEdgeOnLeft)
		{
			double deta = w.disToRoot + leftLen - m_InfoAtVertices[oppositeEdge.indexOfOppositeVert].disUptodate;
			if (deta <= 0)
				return true;
			detaX = xOfVert - w.proportions[0] * edge.length;
			if (detaX * detaX + yOfVert * yOfVert < deta * deta)
				return false;
		}
		else
		{
			double deta = w.disToRoot + rightLen - m_InfoAtVertices[oppositeEdge.indexOfOppositeVert].disUptodate;
			if (deta <= 0)
				return true;
			detaX = xOfVert - w.proportions[1] * edge.length;
			if (detaX * detaX + yOfVert * yOfVert < deta * deta)
				return false;
		}	
	}
	return true;
}
