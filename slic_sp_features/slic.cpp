#include "slic.h"

Dir operator++(Dir &dir, int) {
	Dir result = dir;
	if (dir == Dir::south) {
		dir = Dir::southwest;
	}
	else {
		using IntType = typename underlying_type<Dir>::type;
		dir = static_cast<Dir>(static_cast<IntType>(dir) + 1);
	}
	return result;
}

Dir operator--(Dir &dir, int) {
	Dir result = dir;
	if (dir == Dir::southwest) {
		dir = Dir::south;
	}
	else {
		using IntType = typename underlying_type<Dir>::type;
		dir = static_cast<Dir>(static_cast<IntType>(dir) - 1);
	}
	return result;
}

/*
 * Constructor. Nothing is done here.
 */
Slic::Slic() {

}

/*
 * Destructor. Clear any present data.
 */
Slic::~Slic() {
    clear_data();
}

/*
 * Clear the data as saved by the algorithm.
 *
 * Input : -
 * Output: -
 */
void Slic::clear_data() {
    clusters.clear();
    distances.clear();
    centers.clear();
    center_counts.clear();
}

/*
 * Initialize the cluster centers and initial values of the pixel-wise cluster
 * assignment and distance values.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::init_data(IplImage *image) {
    /* Initialize the cluster and distance matrices. */
    for (int i = 0; i < image->width; i++) { 
        vector<int> cr;
        vector<double> dr;
        for (int j = 0; j < image->height; j++) {
            cr.push_back(-1);
            dr.push_back(FLT_MAX);
        }
        clusters.push_back(cr);
        distances.push_back(dr);
    }
    
    /* Initialize the centers and counters. */
    for (int i = step; i < image->width - step/2; i += step) {
        for (int j = step; j < image->height - step/2; j += step) {
            vector<double> center;
            /* Find the local minimum (gradient-wise). */
            CvPoint nc = find_local_minimum(image, cvPoint(i,j));
            CvScalar colour = cvGet2D(image, nc.y, nc.x);
            
            /* Generate the center vector. */
            center.push_back(colour.val[0]);
            center.push_back(colour.val[1]);
            center.push_back(colour.val[2]);
            center.push_back(nc.x);
            center.push_back(nc.y);
            
            /* Append to vector of centers. */
            centers.push_back(center);
            center_counts.push_back(0);
        }
    }
}

/*
 * Compute the distance between a cluster center and an individual pixel.
 *
 * Input : The cluster index (int), the pixel (CvPoint), and the Lab values of
 *         the pixel (CvScalar).
 * Output: The distance (double).
 */
double Slic::compute_dist(int ci, CvPoint pixel, CvScalar colour) {
    double dc = sqrt(pow(centers[ci][0] - colour.val[0], 2) + pow(centers[ci][1]
            - colour.val[1], 2) + pow(centers[ci][2] - colour.val[2], 2));
    double ds = sqrt(pow(centers[ci][3] - pixel.x, 2) + pow(centers[ci][4] - pixel.y, 2));
    
    return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
    
    //double w = 1.0 / (pow(ns / nc, 2));
    //return sqrt(dc) + sqrt(ds * w);
}

/*
 * Find a local gradient minimum of a pixel in a 3x3 neighbourhood. This
 * method is called upon initialization of the cluster centers.
 *
 * Input : The image (IplImage*) and the pixel center (CvPoint).
 * Output: The local gradient minimum (CvPoint).
 */
CvPoint Slic::find_local_minimum(IplImage *image, CvPoint center) {
    double min_grad = FLT_MAX;
    CvPoint loc_min = cvPoint(center.x, center.y);
    
    for (int i = center.x-1; i < center.x+2; i++) {
        for (int j = center.y-1; j < center.y+2; j++) {
            CvScalar c1 = cvGet2D(image, j+1, i);
            CvScalar c2 = cvGet2D(image, j, i+1);
            CvScalar c3 = cvGet2D(image, j, i);
            /* Convert colour values to grayscale values. */
            double i1 = c1.val[0];
            double i2 = c2.val[0];
            double i3 = c3.val[0];
            /*double i1 = c1.val[0] * 0.11 + c1.val[1] * 0.59 + c1.val[2] * 0.3;
            double i2 = c2.val[0] * 0.11 + c2.val[1] * 0.59 + c2.val[2] * 0.3;
            double i3 = c3.val[0] * 0.11 + c3.val[1] * 0.59 + c3.val[2] * 0.3;*/
            
            /* Compute horizontal and vertical gradients and keep track of the
               minimum. */
            if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3,2)) < min_grad) {
                min_grad = fabs(i1 - i3) + fabs(i2 - i3);
                loc_min.x = i;
                loc_min.y = j;
            }
        }
    }
    
    return loc_min;
}

/*
 * Compute the over-segmentation based on the step-size and relative weighting
 * of the pixel and colour values.
 *
 * Input : The Lab image (IplImage*), the stepsize (int), and the weight (int).
 * Output: -
 */
void Slic::generate_superpixels(IplImage *image, int step, int nc) {
    this->step = step;
    this->nc = nc;
    this->ns = step;
    
    /* Clear previous data (if any), and re-initialize it. */
    clear_data();
    init_data(image);
    
    /* Run EM for 10 iterations (as prescribed by the algorithm). */
    for (int i = 0; i < NR_ITERATIONS; i++) {
        /* Reset distance values. */
        for (int j = 0; j < image->width; j++) {
            for (int k = 0;k < image->height; k++) {
                distances[j][k] = FLT_MAX;
            }
        }

        for (int j = 0; j < (int) centers.size(); j++) {
            /* Only compare to pixels in a 2 x step by 2 x step region. */
            for (int k = centers[j][3] - step; k < centers[j][3] + step; k++) {
                for (int l = centers[j][4] - step; l < centers[j][4] + step; l++) {
                
                    if (k >= 0 && k < image->width && l >= 0 && l < image->height) {
                        CvScalar colour = cvGet2D(image, l, k);
                        double d = compute_dist(j, cvPoint(k,l), colour);
                        
                        /* Update cluster allocation if the cluster minimizes the
                           distance. */
                        if (d < distances[k][l]) {
                            distances[k][l] = d;
                            clusters[k][l] = j;
                        }
                    }
                }
            }
        }
        
        /* Clear the center values. */
        for (int j = 0; j < (int) centers.size(); j++) {
            centers[j][0] = centers[j][1] = centers[j][2] = centers[j][3] = centers[j][4] = 0;
            center_counts[j] = 0;
        }
        
        /* Compute the new cluster centers. */
        for (int j = 0; j < image->width; j++) {
            for (int k = 0; k < image->height; k++) {
                int c_id = clusters[j][k];
                
                if (c_id != -1) {
                    CvScalar colour = cvGet2D(image, k, j);
                    
                    centers[c_id][0] += colour.val[0];
                    centers[c_id][1] += colour.val[1];
                    centers[c_id][2] += colour.val[2];
                    centers[c_id][3] += j;
                    centers[c_id][4] += k;
                    
                    center_counts[c_id] += 1;
                }
            }
        }

        /* Normalize the clusters. */
        for (int j = 0; j < (int) centers.size(); j++) {
            centers[j][0] /= center_counts[j];
            centers[j][1] /= center_counts[j];
            centers[j][2] /= center_counts[j];
            centers[j][3] /= center_counts[j];
            centers[j][4] /= center_counts[j];
        }
    }
}

/*
 * Enforce connectivity of the superpixels. This part is not actively discussed
 * in the paper, but forms an active part of the implementation of the authors
 * of the paper.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::create_connectivity(IplImage *image) {
    int label = 0, adjlabel = 0;
    const int lims = (image->width * image->height) / ((int)centers.size());
    
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    
    /* Initialize the new cluster matrix. */
    vec2di new_clusters;
    for (int i = 0; i < image->width; i++) { 
        vector<int> nc;
        for (int j = 0; j < image->height; j++) {
            nc.push_back(-1);
        }
        new_clusters.push_back(nc);
    }

    for (int i = 0; i < image->width; i++) {
        for (int j = 0; j < image->height; j++) {
            if (new_clusters[i][j] == -1) {
                vector<CvPoint> elements;
                elements.push_back(cvPoint(i, j));
            
                /* Find an adjacent label, for possible use later. */
                for (int k = 0; k < 4; k++) {
                    int x = elements[0].x + dx4[k], y = elements[0].y + dy4[k];
                    
                    if (x >= 0 && x < image->width && y >= 0 && y < image->height) {
                        if (new_clusters[x][y] >= 0) {
                            adjlabel = new_clusters[x][y];
                        }
                    }
                }
                
                int count = 1;
                for (int c = 0; c < count; c++) {
                    for (int k = 0; k < 4; k++) {
                        int x = elements[c].x + dx4[k], y = elements[c].y + dy4[k];
                        
                        if (are_valid_values(image, x, y)) {
                            if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                                elements.push_back(cvPoint(x, y));
                                new_clusters[x][y] = label;
                                count += 1;
                            }
                        }
                    }
                }
                
                /* Use the earlier found adjacent label if a segment size is
                   smaller than a limit. */
                if (count <= lims >> 2) {
                    for (int c = 0; c < count; c++) {
                        new_clusters[elements[c].x][elements[c].y] = adjlabel;
                    }
                    label -= 1;
                }
                label += 1;
            }
        }
    }
}

/*
 * Display the cluster centers.
 *
 * Input : The image to display upon (IplImage*) and the colour (CvScalar).
 * Output: -
 */
void Slic::display_center_grid(IplImage *image, CvScalar colour) {
    for (int i = 0; i < (int) centers.size(); i++) {
        cvCircle(image, cvPoint(centers[i][3], centers[i][4]), 2, colour, 2);
    }
}

/*
 * Display a single pixel wide contour around the clusters.
 *
 * Input : The target image (IplImage*) and contour colour (CvScalar).
 * Output: -
 */
void Slic::display_contours(IplImage *image, CvScalar colour) {
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	/* Initialize the contour vector and the matrix detailing whether a pixel
	 * is already taken to be a contour. */
	vec2db istaken;
	for (int i = 0; i < image->width; i++) { 
        vector<bool> nb;
        for (int j = 0; j < image->height; j++) {
            nb.push_back(false);
        }
        istaken.push_back(nb);
    }
    
    /* Go through all the pixels. */
    for (int i = 0; i < image->width; i++) {
        for (int j = 0; j < image->height; j++) {
            int nr_p = 0;
            
            /* Compare the pixel to its 8 neighbours. */
            for (int k = 0; k < 8; k++) {
                int x = i + dx8[k], y = j + dy8[k];
                
                if (x >= 0 && x < image->width && y >= 0 && y < image->height) {
                    if (istaken[x][y] == false && clusters[i][j] != clusters[x][y]) {
                        nr_p += 1;
                    }
                }
            }
            
            /* Add the pixel to the contour list if desired. */
            if (nr_p >= 2) {
                contours.push_back(cvPoint(i,j));
                istaken[i][j] = true;
            }
        }
    }
    
    /* Draw the contour pixels. */
    for (int i = 0; i < (int)contours.size(); i++) {
        cvSet2D(image, contours[i].y, contours[i].x, colour);
    }
}

/*
 * Give the pixels of each cluster the same colour values. The specified colour
 * is the mean RGB colour per cluster.
 *
 * Input : The target image (IplImage*).
 * Output: -
 */
void Slic::colour_with_cluster_means(IplImage *image) {
    vector<CvScalar> colours(centers.size());
    
    /* Gather the colour values per cluster. */
    for (int i = 0; i < image->width; i++) {
        for (int j = 0; j < image->height; j++) {
            int index = clusters[i][j];
            CvScalar colour = cvGet2D(image, j, i);
            
            colours[index].val[0] += colour.val[0];
            colours[index].val[1] += colour.val[1];
            colours[index].val[2] += colour.val[2];
        }
    }
    
    /* Divide by the number of pixels per cluster to get the mean colour. */
    for (int i = 0; i < (int)colours.size(); i++) {
        colours[i].val[0] /= center_counts[i];
        colours[i].val[1] /= center_counts[i];
        colours[i].val[2] /= center_counts[i];
    }
    
    /* Fill in. */
    for (int i = 0; i < image->width; i++) {
        for (int j = 0; j < image->height; j++) {
            CvScalar ncolour = colours[clusters[i][j]];
            cvSet2D(image, j, i, ncolour);
        }
    }
}

void Slic::colour_superpixels(IplImage *image) {
	vector<CvScalar> colours(centers.size());

	for (int i = 0; i < colours.size(); i++) {
		colours[i] = CvScalar(rand() % 256, rand() % 256, rand() % 256);
	}

	for (int i = 0; i < image->width; i++) {
		for (int j = 0; j < image->height; j++) {
			CvScalar ncolour = colours[clusters[i][j]];
			cvSet2D(image, j, i, ncolour);
		}
	}
}

void Slic::display_vertices(IplImage *image, CvScalar colour) {
	const int dx4[4] = { -1, -1, 0, 0 };
	const int dy4[4] = { 0, -1, -1, 0 }; 

	for (int i = 0; i < image->width; i++) {
		for (int j = 0; j < image->height; j++) {

			if ((j == 0 || j == image->height - 1) && (i == 0 || clusters[i][j] != clusters[i - 1][j]) ||
				(j == 0 || j == image->height - 1) && i == image->width - 1 ||
				(i == 0 || i == image->width - 1) && clusters[i][j] != clusters[i][j - 1]) {
				verts.push_back(cvPoint(i, j));
				continue;
			}

			set<int> adjacent_clusters;

			for (int k = 0; k < 4; k++) {
				int x = i + dx4[k], y = j + dy4[k];

				if (x >= 0 && x < image->width && y >= 0 && y < image->height) {
					adjacent_clusters.insert(clusters[x][y]);
				}
			}

			if (adjacent_clusters.size() >= 3) {
				verts.push_back(cvPoint(i, j));
			}
		}
	}

	for (int i = 0; i < verts.size(); i++) {
		cvSet2D(image, verts[i].y, verts[i].x, colour);
	}
}

void Slic::save_contours(IplImage image, const char* filename) {
	Mat contours_color(image.height, image.width, CV_8UC3, CV_RGB(255, 255, 255)), contours_grayscale, contours_binary;
	for (int i = 0; i < contours.size(); i++) {
		contours_color.at<Vec3b>(contours[i].y, contours[i].x)[0] = 0;
		contours_color.at<Vec3b>(contours[i].y, contours[i].x)[1] = 0;
		contours_color.at<Vec3b>(contours[i].y, contours[i].x)[2] = 0;
	}
	cvtColor(contours_color, contours_grayscale, CV_BGR2GRAY);
	threshold(contours_grayscale, contours_binary, 127, 255, THRESH_BINARY);
	imwrite(filename, contours_binary);
}

void Slic::get_pixel_by_direction(Dir dir, int &new_x, int &new_y, int x, int y) {
	switch (dir) {
	case Dir::southwest:
		new_x = x - 1;
		new_y = y + 1;
		break;
	case Dir::west:
		new_x = x - 1;
		new_y = y;
		break;
	case Dir::northwest:
		new_x = x - 1;
		new_y = y - 1;
		break;
	case Dir::north:
		new_x = x;
		new_y = y - 1;
		break;
	case Dir::northeast:
		new_x = x + 1;
		new_y = y - 1;
		break;
	case Dir::east:
		new_x = x + 1;
		new_y = y;
		break;
	case Dir::southeast:
		new_x = x + 1;
		new_y = y + 1;
		break;
	case Dir::south:
		new_x = x;
		new_y = y + 1;
		break;
	}
}

bool Slic::are_valid_values(IplImage *image, int new_x, int new_y) {
	if (new_x >= 0 && new_x < image->width && new_y >= 0 && new_y < image->height) return true;
	return false;
}

bool Slic::is_rotation_needed(IplImage *image, int new_x, int new_y, int x, int y) {
	if (!are_valid_values(image, new_x, new_y)) {
		return true;
	}

	if (is_vertex[new_x][new_y]) {
		return false;
	}

	if (is_edge[new_x][new_y] >= 2) {
		cout << "!!" << endl;
		return true;
	}

	if (clusters[new_x][new_y] != clusters[x][y]) {
		return true;
	}

	return false;
}

void Slic::go_next_edge_pixel(IplImage *image, Dir &dir, int &x, int &y) {
	Dir init_dir = dir;
	int new_x, new_y;
	get_pixel_by_direction(dir, new_x, new_y, x, y);
	while (is_rotation_needed(image, new_x, new_y, x, y)) {
		/*dir = rotate_ccw(dir);*/
		dir--;
		if (dir == init_dir) { 
			return;
		}
		get_pixel_by_direction(dir, new_x, new_y, x, y);
	}
	x = new_x, y = new_y;
}

void Slic::construct_graph(IplImage *image) {

	cout << 1 << endl; 

    typedef pair<CvPoint, CvPoint> Edge;

	for (int i = 0; i < image->width; i++) {
		vector<bool> column;
		for (int j = 0; j < image->height; j++) {
			column.push_back(false);
		}
		is_vertex.push_back(column);
	}

	for (int i = 0; i < verts.size(); i++) {
		is_vertex[verts[i].x][verts[i].y] = true;
	}

	for (int i = 0; i < image->width; i++) {
		vector<int> column;
		for (int j = 0; j < image->height; j++) {
			column.push_back(0);
		}
		is_edge.push_back(column);
	}

	vector<Edge> edges;
	vector<int> weights;

	for (int i = 0; i < verts.size(); i++) {
		int x = verts[i].x, y = verts[i].y;
		int weight_counter = 0;

		Dir dir = Dir::southwest;
		while (true) {
			dir++;
			go_next_edge_pixel(image, dir, x, y);
			if (is_vertex[x][y]) { break; }
			is_edge[x][y]++;
			cvSet2D(image, y, x, CV_RGB(255, 255, 0));
			weight_counter++;
		}
		
		if (x == verts[i].x && y == verts[i].y) { continue; }

		edges.push_back(Edge(verts[i], CvPoint(x, y)));
		weights.push_back(weight_counter);
	}
}
