#ifndef SLIC_H
#define SLIC_H

/* slic.h.
 *
 * Written by: Pascal Mettes.
 *
 * This file contains the class elements of the class Slic. This class is an
 * implementation of the SLIC Superpixel algorithm by Achanta et al. [PAMI'12,
 * vol. 34, num. 11, pp. 2274-2282].
 *
 * This implementation is created for the specific purpose of creating
 * over-segmentations in an OpenCV-based environment.
 */

#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <set>
#include <map>

#include <iostream>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
using namespace boost;

using namespace std;
using namespace cv;

/* 2d matrices are handled by 2d vectors. */
#define vec2dd vector<vector<double> >
#define vec2di vector<vector<int> >
#define vec2db vector<vector<bool> >
/* The number of iterations run by the clustering algorithm. */
#define NR_ITERATIONS 10

enum class Dir { 
	southwest = 0, 
	west = 1,
	northwest = 2,
	north = 3,
	northeast = 4,
	east = 5,
	southeast = 6,
	south = 7	
};

/*
 * class Slic.
 *
 * In this class, an over-segmentation is created of an image, provided by the
 * step-size (distance between initial cluster locations) and the colour
 * distance parameter.
 */
class Slic {
    private:
        /* The cluster assignments and distance values for each pixel. */
        vec2di clusters;
        vec2dd distances;
        
        /* The LAB and xy values of the centers. */
        vec2dd centers;
        /* The number of occurences of each center. */
        vector<int> center_counts;
        
        /* The step size per cluster, and the colour (nc) and distance (ns)
         * parameters. */
        int step, nc, ns;

		vector<CvPoint> contours;
		vector<CvPoint> verts;

		vec2db is_vertex;
		vec2di is_edge;

		enum edge_t {left, right};

		typedef pair<CvPoint, CvPoint> Edge;
		typedef pair<int, int> Edge1;

		vector<Edge> edges;
		vector<Edge1> edges1;
		vector<int> weights;
        
        /* Compute the distance between a center and an individual pixel. */
        double compute_dist(int ci, CvPoint pixel, CvScalar colour);
        /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
        CvPoint find_local_minimum(IplImage *image, CvPoint center);
        
        /* Remove and initialize the 2d vectors. */
        void clear_data();
        void init_data(IplImage *image);

		void go_next_edge_pixel(edge_t edge_type, IplImage *image, Dir &dir, int &x, int &y, int weight);
		void get_pixel_by_direction(Dir dir, int &new_x, int &new_y, int x, int y);
		bool are_valid_values(IplImage *image, int new_x, int new_y);
		bool is_rotation_needed(IplImage *image, int new_x, int new_y, int x, int y, int weight);
		void find_edge(edge_t edge_type, IplImage *image, int i);

    public:
        /* Class constructors and deconstructors. */
        Slic();
        ~Slic();
        
        /* Generate an over-segmentation for an image. */
        void generate_superpixels(IplImage *image, int step, int nc);
        /* Enforce connectivity for an image. */
        void create_connectivity(IplImage *image);
        
        /* Draw functions. Resp. displayal of the centers and the contours. */
        void display_center_grid(IplImage *image, CvScalar colour);
        void display_contours(IplImage *image, CvScalar colour);
        void colour_with_cluster_means(IplImage *image);
		void colour_superpixels(IplImage *image);
		void display_vertices(IplImage *image, CvScalar colour);

		void save_contours(IplImage image, const char* filename);
		void construct_graph(IplImage *image);
};

#endif
