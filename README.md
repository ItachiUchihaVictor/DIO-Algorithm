# DIO-Algorithm

The project provides the implementation of the  “On Efficient Shortest Path Computation on Terrain Surface: A Direction-Oriented Approach” stated in Reference. Please kindly cite the paper in the Reference if you used our code. 

# Reference

“On Efficient Shortest Path Computation on Terrain Surface: A Direction-Oriented Approach”, Victor Junqiu Wei, Raymond Chi-Wing Wong, Cheng Long, David M. Mount, Hanan Samet, IEEE Transaction on Knowledge and Data Engineering (TKDE), 2024

# How to Compile and Run the Code 

g++ -o main main_dist.cpp -std=c++11

./main small_terrain.off

Our DIO algorithm is implemented in the function shortestpath_GB(). 

# How to Prepare the Datasets

First download the file through the following link:

[https://www.dropbox.com/s/ofa9ddk138x91w3/dataset.tar.gz?dl=0](https://www.dropbox.com/scl/fi/ykutpuruep1996i1gx51k/TerrainDataset.zip?rlkey=9mm5b0qy4o59s5isx7vj49jts&dl=0)

Data Format:

We used the .off format in the experiment. The content of the .off file is as follows: 

OFF

Number_of_vertices Number_of_faces Number_of_edges

x_coordinate_of_1st_vertex y_coordinate_of_1st_vertex z_coordinate_of_1st_vertex

x_coordinate_of_2nd_vertex y_coordinate_of_2nd_vertex z_coordinate_of_2nd_vertex

......

x_coordinate_of_the_last_vertex y_coordinate_of_the_last_vertex z_coordinate_of_the_last_vertex

ID_of_the_1st_vertex_of_the_1st_face ID_of_the_2nd_vertex_of_the_1st_face ID_of_the_3td_vertex_of_the_1st_face

ID_of_the_1st_vertex_of_the_2nd_face ID_of_the_2nd_vertex_of_the_2nd_face ID_of_the_3td_vertex_of_the_2nd_face

......

ID_of_the_1st_vertex_of_the_last_face ID_of_the_2nd_vertex_of_the_last_face ID_of_the_3td_vertex_of_the_last_face

Each .off data could be visualized by the terrain tool (http://rwcpu1.cse.ust.hk/terrain/).
