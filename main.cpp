#include "distance.h"
#include<sstream>
#include<unistd.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <iostream>
#include <sstream>
#include "geodesic_algorithm_exact.h"
#define qtimes 100
#define lsize 1 
//-------------------------------------
// MAIN
//-------------------------------------
int query_type=0;
int clpsize=6;
int algo_type=0;
int k=3;
int i;
int x[upper_poi];
FILE *fp;
char prefix[255];
std::vector<geodesic::GeodesicAlgorithmExact*> landmarks;
//geodesic::GeodesicAlgorithmExact algorithm;	//create exact algorithm for the mesh
#ifndef WIN32
    double Time_preprocess=0;
    double  Time_dquery=0, Time_knnquery, Time_clpquery;
    double  Space_preprocess=0;
    double  Space_query=0;
    double errorbound_dis, errorbound_knn=0;
    struct rusage myTime_program_start, myTime_preprocess_end, myTime_query_begin, myTime_query_end;
#endif

double shortestpath_GB(int x, int y, geodesic::GeodesicAlgorithmExact& shortestpath){
    geodesic::SurfacePoint source(&mesh.vertices()[x]);
    geodesic::SurfacePoint dest(&mesh.vertices()[y]);
    std::vector<geodesic::SurfacePoint> sources;
    std::vector<geodesic::SurfacePoint> dests;
    sources.clear();
    dests.clear();
    sources.push_back(source);
    dests.push_back(dest);
    shortestpath.propagate_GB(sources, &dests);
    double dist;
    shortestpath.best_source(dest, dist);
    return dist;
}
double shortestpath_LA(int x, int y, geodesic::GeodesicAlgorithmExact& shortestpath, std::vector<geodesic::GeodesicAlgorithmExact *> * landmarks){
    geodesic::SurfacePoint source(&mesh.vertices()[x]);
    geodesic::SurfacePoint dest(&mesh.vertices()[y]);
    std::vector<geodesic::SurfacePoint> sources;
    std::vector<geodesic::SurfacePoint> dests;
    sources.clear();
    dests.clear();
    sources.push_back(source);
    dests.push_back(dest);
    shortestpath.propagate_LA(sources, landmarks,  &dests);
    double dist;
    shortestpath.best_source(dest, dist);
    return dist;
}
double shortestpath_MMP(int x, int y, geodesic::GeodesicAlgorithmExact& shortestpath){
    geodesic::SurfacePoint source(&mesh.vertices()[x]);
    geodesic::SurfacePoint dest(&mesh.vertices()[y]);
    std::vector<geodesic::SurfacePoint> sources;
    std::vector<geodesic::SurfacePoint> dests;
    sources.clear();
    dests.clear();
    sources.push_back(source);
    dests.push_back(dest);
    shortestpath.propagate(sources, 0.0,  &dests);
    double dist;
    shortestpath.best_source(dest, dist);
    return dist;
}
double shortestpath_MMP(int x, geodesic::GeodesicAlgorithmExact& shortestpath, int y){
    geodesic::SurfacePoint source(&mesh.vertices()[x]);
    geodesic::SurfacePoint dest(&mesh.vertices()[y]);
    std::vector<geodesic::SurfacePoint> sources;
    sources.clear();
    sources.push_back(source);
    shortestpath.propagate(sources, geodesic::GEODESIC_INF);
    double dist;
    shortestpath.best_source(dest, dist);
    return dist;
}

int main(int argc, char **argv) 
{
	if(argc < 2)
	{
		std::cout << "usage: mesh_file_name " << std::endl; //try: "hedgehog_mesh.txt 3 14" or "flat_triangular_mesh.txt 1"
		return 0;
	}

 //   s = atof(argv[2]);
	bool success = geodesic::read_mesh_from_file(argv[1],points,faces);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return 0;
	}

    strcpy(prefix, argv[1]);
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh

    landmarks.resize(lsize);
    int x[lsize];

    for(int i=0;i<lsize;i++){
        x[i] = rand()*rand()%mesh.vertices().size();
        for(int j=0;j<i;j++){
            if(i==0)break;
            if(x[i]==x[j]){
                i--;
                break;
            }
        }
    }

    for(int i=0;i<lsize;i++){
        landmarks[i] = new geodesic::GeodesicAlgorithmExact(&mesh);
        geodesic::SurfacePoint p(&mesh.vertices()[x[i]]);
        std::vector<geodesic::SurfacePoint> sources;
        sources.clear();
        sources.push_back(p);
        landmarks[i]->propagate(sources);
    }

    double begin, end;
   double GBtime, LAtime, MMPtime, MMPstime;
   double GBt, LAt, MMPt, MMPst, LAinterval, LAedges;
   double GBinterval, MMPsinterval, GBedges, GBtotaledges;
   double GBtotalint, MMPstotalint, MMPedges, MMPtotaledges;
   double distance;

//   std::ofstream output(std::string(argv[1])+"_Sp.txt", std::ios::out );
   std::ofstream GB(std::string(argv[1])+"_gb.txt", std::ios::out );
//   std::ofstream LA(std::string(argv[1])+"_la.txt", std::ios::out | std::ios::app);
//   std::ofstream MMP(std::string(argv[1])+"_mmp.txt", std::ios::out | std::ios::app);
//   std::ofstream MMPs(std::string(argv[1])+"_mmps.txt", std::ios::out );
   
    for(int i=0;i<qtimes;i++){
        int src = rand()*rand()%mesh.vertices().size();
        int dst = rand()*rand()%mesh.vertices().size();
        
        while(src==dst){
            dst = rand()*rand()%mesh.vertices().size();
        }
        begin = clock(); 
        distance = shortestpath_GB(src, dst, algorithm);
        end = clock();
        GBtime = (end - begin)*1000.0/CLOCKS_PER_SEC;
        algorithm.print_statistics(GBinterval, GBedges);
        GB << distance << " " << GBtime << " " << GBinterval << " " << GBedges/2.0 << std::endl;
        GBt+=GBtime;
        GBtotalint+=GBinterval;
        GBtotaledges+=GBedges/2.0;
/*        begin = clock();
        distance = shortestpath_LA(src, dst, algorithm,  &landmarks);
        end = clock();
        LAtime = (end - begin)*1000.0/CLOCKS_PER_SEC;
        algorithm.print_statistics(LAinterval, LAedges);
        LA << distance << " " << LAtime << " " << LAinterval << std::endl;
        LAt+=LAtime;
*/
/*
        begin = clock();
        //shortestpath_MMP(1000%mesh.vertices().size(), 100000%mesh.vertices().size(), algorithm);
        distance = shortestpath_MMP(src, algorithm, dst);
        end = clock();
        MMPtime = (end - begin)*1000.0/CLOCKS_PER_SEC;
        MMP << distance << " " << MMPtime << std::endl;
        MMPt+=MMPtime;
        begin = clock();
        //shortestpath_MMP(1000%mesh.vertices().size(), 100000%mesh.vertices().size(), algorithm);
        distance = shortestpath_MMP(src, dst, algorithm);
        end = clock();
        MMPstime = (end - begin)*1000.0/CLOCKS_PER_SEC;
        algorithm.print_statistics(MMPsinterval,MMPedges );
        MMPs << distance << " " << MMPstime << " " << MMPsinterval << " " << MMPedges/2.0 << std::endl;
        MMPst+=MMPstime;
        MMPstotalint+=MMPsinterval;
        MMPtotaledges+=MMPedges/2.0;
*/
    }
        
        GBt/=qtimes;
        GBtotalint/=qtimes;
        GBtotaledges/=qtimes;
 //       LAt/=qtimes;
        MMPt/=qtimes;
        MMPst/=qtimes;
        MMPstotalint/=qtimes;
        MMPtotaledges/=qtimes;

    
//        output << GBt << " " << LAt << " " << MMPt << " " << MMPst << std::endl;
 //       output << GBt << " " << GBtotalint << " " << GBtotaledges << " " << MMPst << " " << MMPstotalint << " " << MMPtotaledges << std::endl;
        GB.close();
//        LA.close();
//        MMPs.close();
//        output.close();


//    std::cout << s << std::endl;

	//geodesic::GeodesicAlgorithmSubdivision subdivision_algorithm(&mesh,2);	//with subdivision_level=0 this algorithm becomes Dijkstra, with subdivision_level->infinity it becomes exact
    // WRITE THE RANDOM POINT FILE 
/*    fp = fopen("POINT.C","w");
    if ( fp == NULL )
    {
        puts ( "Cannot open file" );
        exit(1);
    }*/
//    x[0]=randn(mesh.vertices().size());
//    for(i=1;i<poi;i++)
//    {
      //  bool key=true;
     //   while(key){
        //    key=false;
  //          x[i]=(x[i-1]+1)%(mesh.vertices().size());
     /*       for(int j=0;j<i;j++){
                if(x[i]==x[j]){
                    key=true;
                    break;
                }
            }*/
      //  }
 //       fprintf(fp,"%d\n",x[i]);
   // }
 //   fclose(fp);
// READ THE RANDOM POINT FILE AND ASSIGN TO ROOT Node
/*    std::cout<<"Tree Building:"<<std::endl;
    std::cout<<"1: 2D QuadTree"<<std::endl;
    std::cout<<"2: Geodesic Tree"<<std::endl;
    std::cout<<"3: Graph Tree"<<std::endl;
    fp=fopen("output.txt","w");
    fclose(fp);
    for(algo_type=4;algo_type<5;algo_type++){
    pairs=0;
    quadpairs.clear();
    pairvector.clear();
    nodevector.clear();
    geopairs.clear();
    geopairsvector.clear();
    geonodevector.clear();
    graphpairs.clear();
    graphpairsvector.clear();
    graphnodevector.clear();
    std::cout<<"----------Algorithm "<<algo_type<<"------------"<<std::endl;

    if(algo_type==1||algo_type==2)algo_1_2(algorithm);
    if(algo_type==3){
        algo_3(algorithm);
    }
    if(algo_type==4){
        algo_4(algorithm);
    }
    }*/
}
