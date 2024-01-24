#include "distance.h"
#include<cmath>
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

    double MaxDistance = 0;
    double MinDistance = 1000000000;
    double radius=0;
    double distance;
    double const distance_limit = 0;
    geodesic::GeodesicAlgorithmExact landmark(&mesh);
    geodesic::SurfacePoint p(&mesh.vertices()[0]);
    geodesic::SurfacePoint source(&mesh.vertices()[0]);
    geodesic::SurfacePoint target(&mesh.vertices()[1]);
    std::vector<geodesic::SurfacePoint> stop_points(1, target);
    std::vector<geodesic::SurfacePoint> all_sources(1,source);

	for(unsigned i=1; i<mesh.vertices().size(); ++i){
            stop_points.push_back(geodesic::SurfacePoint(&mesh.vertices()[i]));
        }
	    landmark.propagate(all_sources, distance_limit, &stop_points);
	for(unsigned i=1; i<mesh.vertices().size(); ++i){
            geodesic::SurfacePoint p(&mesh.vertices()[i]);
            
            landmark.best_source(p,distance);
            radius= std::max(distance, radius);
            MinDistance= std::min(distance, MinDistance);
            std::cout << "radius: " << distance << std::endl;
        }
    MaxDistance=radius; 
//    MaxDistance*=2; 
    std::cout << "Max distance obtained. " << std::endl;
    std::cout << "max distance is " << MaxDistance << std::endl;
    std::cout << "Min distance obtained. " << std::endl;
    std::cout << "min distance is " << MinDistance << std::endl;

    double begin, end;
   double GBtime, LAtime, MMPtime, MMPstime;
   double GBt, LAt, MMPt, MMPst, LAinterval, LAedges;
   double GBinterval, MMPsinterval, GBedges, GBtotaledges;
   double GBtotalint, MMPstotalint, MMPedges, MMPtotaledges;

//   std::ofstream output(std::string(argv[1])+"_Sp.txt", std::ios::out );
   std::vector<std::vector<int>> sources;
   std::vector<std::vector<int>> destinations;
   sources.resize(10);
   destinations.resize(10);
//   std::ofstream LA(std::string(argv[1])+"_la.txt", std::ios::out | std::ios::app);
//   std::ofstream MMP(std::string(argv[1])+"_mmp.txt", std::ios::out | std::ios::app);
//   std::ofstream MMPs(std::string(argv[1])+"_mmps.txt", std::ios::out );
   
    while(true){
        bool complete = true;
        for(i=0;i<10;i++){
            if(sources[i].size() < 100){ complete = false;}
            if(destinations[i].size() < 100){ complete = false;}
            std::cout << "source[" << i << "] size: " << sources[i].size() << " destinations[" << i << "].size(): " << destinations[i].size() << std::endl;
        }
        if(complete){ break; }
        int src = rand()*rand()%mesh.vertices().size();
        
        geodesic::SurfacePoint source(&mesh.vertices()[src]);
        geodesic::SurfacePoint target(&mesh.vertices()[1]);
        std::vector<geodesic::SurfacePoint> stop_points(1, target);
        std::vector<geodesic::SurfacePoint> all_sources(1,source);

	for(unsigned i=0; i<mesh.vertices().size(); ++i){
            stop_points.push_back(geodesic::SurfacePoint(&mesh.vertices()[i]));
        }
	    landmark.propagate(all_sources, distance_limit, &stop_points);
        int numGenPoints = 0;
	for(unsigned i=0; i<mesh.vertices().size(); ++i){
            int dst = i; //rand()*rand()%mesh.vertices().size();
  //      while(src==dst){
  //          dst = rand()*rand()%mesh.vertices().size();
  //      }
            if(numGenPoints>100) break;
            if(dst==src){continue;}
            geodesic::SurfacePoint p(&mesh.vertices()[dst]);
            
            landmark.best_source(p,distance);
        std::cout << "distance: " << distance << " MinDistance: " << MinDistance << " MaxDistance: " << MaxDistance << std::endl;
        //int j = log10(distance/MinDistance)/log10(2);
        int j = (distance-MinDistance)/((MaxDistance-MinDistance)/10.0);
        j = std::min(j, 9);
         
        std::cout << "j: " << j << " src: " << src << " dst: " << dst << std::endl;
        if(sources[j].size() < 100){
            sources[j].push_back(src);
            numGenPoints++;}
        if(destinations[j].size() < 100)destinations[j].push_back(dst);
        }
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
    std::ofstream GB(std::string(argv[1])+"_query.txt", std::ios::out );
    for(int i=0;i<10;i++){
        GB << i << std::endl;
        for(int j=0;j<sources[i].size();j++) GB << sources[i][j] << " " << destinations[i][j] << std::endl;
    }    
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
