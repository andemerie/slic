/*
 * test_slic.cpp.
 *
 * Written by: Pascal Mettes.
 *
 * This file creates an over-segmentation of a provided image based on the SLIC
 * superpixel algorithm, as implemented in slic.h and slic.cpp.
 */

#include "slic.h"
#include <fstream>

ofstream fout("report.txt");

int main(int argc, char* argv[]) {
  string images[] = {
      "dog.png", "lena.bmp", "forest.jpg", "earth.jpg", "doge.jpg"
  };

  for (string image_name : images) {
    /* Load the image and convert to Lab colour space. */
    string image_path = "test/" + image_name;
    IplImage* image = cvLoadImage(image_path.c_str(), 1);
    IplImage* lab_image = cvCloneImage(image);
    cvCvtColor(image, lab_image, CV_BGR2Lab);

    /* Yield the number of superpixels and weight-factors from the user. */
    int w = image->width, h = image->height;
    int nr_superpixels = 400;
    int nc = 40;

    double step = sqrt((w * h) / (double) nr_superpixels);

    /* Perform the SLIC superpixel algorithm. */
    Slic slic;
    auto start_time = std::chrono::steady_clock::now();
    slic.generate_superpixels(lab_image, step, nc);
    slic.create_connectivity(lab_image);

    /* Display the contours and show the result. */
    slic.display_contours(image, CV_RGB(255, 0, 0));
//    cvShowImage("result", image);
//    cvWaitKey(0);

    auto end_time = std::chrono::steady_clock::now();
    fout << "image: " << image_name << endl;
    fout << "size: " << w << "x" << h << " = " << w * h << " pixels" << endl;
    fout << "elapsed time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(
             end_time - start_time
         ).count() / 1000. << "s" << endl;;
    fout << "---------------------------------------------" << endl << endl;


    string out_path = "out/" + image_name;
    cvSaveImage(out_path.c_str(), image);
  }
}
