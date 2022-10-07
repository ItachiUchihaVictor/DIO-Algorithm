/* 
 * File:   geodesic.h
 * Author: mkaul
 *
 * Created on March 12, 2014, 8:48 AM
 */
#include <vector>
#include <algorithm>
#include <list>
#include <queue>
#include <map>
#include <iterator>

#include <iostream>
#include <sstream>
// #include <chrono>

#include <string>
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <cstring>
#include <exception>
#include <stdexcept>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>                                            /* For command line args */
#include <boost/config.hpp>
#include <boost/algorithm/string.hpp>


#ifndef GEODESIC_H
#define	GEODESIC_H
#define TIMING



using namespace std;
//using namespace boost;

class Edge; 

/* 
 * Vertex/Node 
 */
class Node 
{
public:
	int           id;
  double        x;
  double        y;
  double        z;
	list<Edge*>   neighbors;
  vector<int>   adj_face_ids;                                                   /* Needed to know the faces that share an edge */
	bool          visited;
	double        cost;
	Node          *back_pointer;

	Node(int id): id(id), visited(false), cost(0), x(0), y(0), z(0){}
};

/* 
 * Edge 
 */
class Edge 
{
public:
	double weight;
	Node   *dest;

	Edge(double c, Node *d = NULL): weight(c),dest(d){}
	~Edge(){if(dest)delete dest;}
};

/* 
 * Node Map 
 */
class NodeMap{
public:
	map<int,Node*> nodemap;

	Node* find_in_nodemap(int id){
		Node *result = nodemap[id];
		if ( result == 0 )
			result = nodemap[id] = new Node(id);
		return result;
	}
  
	friend ostream& operator<<(ostream& o, NodeMap nm){
	
    map<int,Node*>::iterator im;
		for(im = nm.nodemap.begin(); im != nm.nodemap.end(); im++){
			o << "*****  Node: " << (*im).second->id 
                    << "("  << (*im).second->x <<","
                            << (*im).second->y <<","
                            << (*im).second->z <<")"  
                    << endl;
	    list<Edge*> neighbors = (*im).second->neighbors;
			list<Edge*>::iterator e;
			for(e = neighbors.begin(); e != neighbors.end(); e++){
				cout << "   -> " << (*e)->dest->id << " weight " << (*e)->weight <<endl;
			}
		}
		return o;
	} /* End of friend */
  
};

/* 
 * Edge De-duplication Map 
 */
class EdgeMap{
public:
	map< std::pair<int,int>, int> edgemap;

	double find_in_edgemap(int a, int b){
	  return edgemap[ std::make_pair(a,b) ];
	}
  
  void add_to_edgemap(int a, int b, double w){
    //cout << "add_to_edge_map: "<< a << "," << b << ")\n";
    edgemap[ std::make_pair(a,b) ] = w;
    edgemap[ std::make_pair(b,a) ] = w;
  }
  
  friend ostream& operator<<(ostream& o, EdgeMap em){
	
    map< std::pair<int,int>, int>::iterator im;
		for(im = em.edgemap.begin(); im != em.edgemap.end(); im++){
       o << "*****  Edge: ("  << (*im).first.first << "," 
                              << (*im).first.second <<") --> "
                              << (*im).second << "\n";  
    }
    return o;
  }
  
};


/*
#ifdef TIMING
#define INIT_TIMER auto st = std::chrono::high_resolution_clock::now();
#define START_TIMER     st = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-st \
    ).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif
*/



#endif	/* GEODESIC_H */

