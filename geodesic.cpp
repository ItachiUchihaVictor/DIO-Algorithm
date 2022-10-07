// geodesic.cpp : 
//

#include <cstdio>
#include <iostream>
#include <iomanip>
#include<fstream>
#include<sstream>
#include<unistd.h>

#include <cstdlib>
#include <assert.h>

#include <algorithm>
#include <vector>
#include <map>
#include <iterator>

#include<cmath>
#include <sys/time.h>

#include "BaseModel.h"
#include "RichModel.h"
#include "ExactMethodForDGP.h"
#include "PreviousCH.h"
 //#include "ICHWithFurtherPriorityQueue.h"
#include "ImprovedCHWithEdgeValve.h"
#include "geodesic.h"

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>   //  For command line args 
#include <boost/config.hpp>

//#include <Windows.h>
#define qtimes 500
using namespace std;

double GetTime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}



int main(int argc, char* argv[])
{
	std::string inf;
        std::string outf;
  int S=0,T=0;
  NodeMap nodes;
  EdgeMap edges;
  
    std::ofstream resultfile(std::string(argv[1])+"_dist_AstarICHsp.txt", std::ios::out);
    std::ofstream distICH(std::string(argv[1])+"_dist_AstarICH.txt", std::ios::out);
//    std::ofstream distCH(std::string(argv[1])+"_CH.txt", std::ios::out);
//  resultfile.open(outf);
 /* 
  resultfile << "----------------------Algorithm begins:----------------\n";
        resultfile << "Algorithm name:\t" << alg.GetAlgorithmName() << endl;
	alg.Execute();
  
  resultfile << "Running time (secs):\t" << alg.GetRunTime()/CLOCKS_PER_SEC <<  endl;
	resultfile << "Peak memory: (M)\t"  << alg.GetMemoryCost() << endl;
*/
    double CHt, ICHt, CHseg, CHver, ICHseg, ICHver;
    double totalCH=0, totalICHt=0, totalCHseg=0, totalCHver=0, totalICHseg=0, totalICHver=0;
    double begin, end;

  // Declare the supported options.
  /*
  boost::program_options::options_description d("Allowed options for geodesic.cpp");
  d.add_options()
      ("help","produce this help message")
      ("off",  boost::program_options::value<std::string>(), "load graph from file 'arg'")
      ("start",boost::program_options::value<int>(),         "start = 'arg'")
      ("end",boost::program_options::value<int>(),           "end = 'arg'")
      ("output", boost::program_options::value<std::string>(), "write result to file 'arg'")
      ;

  boost::program_options::variables_map m;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, d), m);
  boost::program_options::notify(m);
*/
/*  if (m.count("help")){
    //Display the options_description
    std::cout << d << "\n";
    std::cout << "./geodesic --in something.off --start 0 --end 100 \n";
  }
  if(m.count("output")){
      outf = m["output"].as<std::string>();
      std::cout << "output result file: '" << outf << "'\n";
  }
  if(m.count("off") ){
    inf = m["off"].as<std::string>();
    std::cout << "input off file : '"<< inf << "'\n";
  }
  if(m.count("start") ){
    S = m["start"].as<int>();
    std::cout << "start vertex : '"<< S << "'\n";
  }
  if(m.count("end") ){
    T = m["end"].as<int>() ;
    std::cout << "destination vertex : '"<< T << "'\n";
  }
  else{
    std::cout << "No args set\n";
  }
  */
  
  /* JUST TESTING !!! */
    /* Reload the model */
//  CRichModel remodel(inf.c_str());
        CRichModel remodel(argv[1]);
	cout <<"----------------------model info begin----------------\n";	
	cout << "Tmp File name:\t" << remodel.GetFileName() << endl;		
	try
	{
		remodel.LoadModel();
	}
	catch(const char* msg)
	{
		cout << "ERRORS happen!\n" << msg << endl;		
		return 1;
	}
  
	remodel.Preprocess();
  cout <<"\n\n Model was preprocessed too\n"; 

	if (!remodel.HasBeenLoad() || !remodel.HasBeenProcessed())
	{
		cout << "The model fails to be handled." << endl;
		return 1;
	}
  
  cout << "Face number:\t" << remodel.GetNumOfFaces() << endl;
	cout << "Vertex number:\t" << remodel.GetNumOfVerts() << endl;
	cout << "Edge number:\t" << remodel.GetNumOfEdges() << endl;
	if (remodel.IsClosedModel())
		cout << "It is a closed model." << endl;
	else 
		cout << "It is an open model with holes." << endl;
	if (remodel.GetNumOfIsolated() > 1)
	{
		cout << "Perhaps it is composed of several components." << endl;
	}
	cout <<"----------------------model info end------------------\n";

  
     std::ofstream querypair("query.txt", std::ios::out); 
	set<int> sources;
	// sources.insert(0);//the vertex id is zero based.
	//here you can add more source points...
	//sources.insert(1);
//	sources.insert(0);
	// CICHWithFurtherPriorityQueue alg(remodel, sources);
   int srcs[10][100];
   int dsts[10][100];
   std::ifstream query(std::string(argv[1])+"_query.txt", std::ios::out );
   for(int i=0;i<10;i++){
       int num;
       query >> num;
 //      std::cout << "query " << num << ":" << std::endl;  
       for(int j=0;j<100;j++){
           query >> srcs[i][j] >> dsts[i][j];
       }
   }
//   for(int i=0;i<10;i++){
//       for(int j=0;j<100;j++){
//           std::cout << srcs[i][j] << " " << dsts[i][j] << std::endl;
//       }
//   }
    for(int i=0;i<10;i++){
      totalICHt = 0;
      totalICHseg = 0;
      totalICHver = 0;
      for(int j=0;j<100;j++){
//        std::cout << "src: " << srcs[i][j] << " dst: " << dsts[i][j] << std::endl; 
        S = srcs[i][j]; //std::abs(rand()*rand()%remodel.GetNumOfVerts());
        sources.insert(S);
        T = dsts[i][j]; // std::abs(rand()*rand()%remodel.GetNumOfVerts());
        querypair << S << " " << T << std::endl;
   //     std::cout << "src: " << S << " dst: " << T << std::endl; 
//        while(T==S)T = std::abs(rand()*rand()%remodel.GetNumOfVerts());
     //   std::cout << "S T: " << S << " " << T << std::endl;
        CImprovedCHWithEdgeValve ICHalg(remodel, sources);
        CPreviousCH CHalg(remodel, sources);
        ICHalg.SetDest(T);

        begin = GetTime(); //clock();
        ICHalg.Execute();
        end = GetTime(); //clock();
        ICHt = (-begin+end)*1000.0;
        totalICHt+=ICHt;
       // if(ICHalg.GetDistanceAt(T)<6000){
//		i--;
//	 	continue;	
//	}
        distICH << ICHalg.GetDistanceAt(T) << " " << ICHt << " " << ICHalg.GetWindowNum() << " " << ICHalg.visited_vertices << std::endl;
        totalICHseg += ICHalg.GetWindowNum(); //GetMaxLenOfQue();
        totalICHver += ICHalg.visited_vertices;
/*
        begin = clock();
        CHalg.Execute();
        end =clock();
        CHt = (-begin+end)*1000.0/CLOCKS_PER_SEC;
        totalCH+=CHt;
        distCH << CHalg.GetDistanceAt(T) << " " << CHt << " " << CHalg.GetMaxLenOfQue() << " " << remodel.GetNumOfVerts() << std::endl;
        totalCHseg += CHalg.GetMaxLenOfQue();
*/
    //    delete &CHalg;
    //    delete &ICHalg;
        sources.clear();
    }
 
   resultfile << i << " " 
       <<  totalICHt/qtimes <<  " " << totalICHseg/qtimes << " " << totalICHver/qtimes << " "
       //<<  totalCH/qtimes <<  " " << totalCHseg/qtimes << " " << remodel.GetNumOfVerts() 
       << std::endl;
  }
    resultfile.close();
//    distCH.close();
    distICH.close();

  //alg.Execute();
	//The following line gets the geodesic distance at the vertex 10 
//	double dis = alg.GetDistanceAt(100);
	//
//	vector<EdgePoint> path;
	//sometimes when the source ids are not unique, the following 
	// code can return the source that is nearest to the destination.
	// at the same, the code gets the path
//	int source = alg.FindSourceVertex(100, path);
	
  //double x=0,y=0,z=0;
/*  
  resultfile << dis << " " << path.size() << endl;
  
  for(int i=0; i < path.size(); i++){
    resultfile  << path[i].Get3DPoint(remodel).x << " " 
                << path[i].Get3DPoint(remodel).y << " "
                << path[i].Get3DPoint(remodel).z << endl;
  }
  */
  return 0;
  
 
}

