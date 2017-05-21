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
    auto start_time = std::chrono::steady_clock::now();
    /* Load the image and convert to Lab colour space. */
    Mat image = imread("test/" + image_name, 1);
    Mat lab_image;
    cvtColor(image, lab_image, COLOR_BGR2Lab);

    /* Yield the number of superpixels and weight-factors from the user. */
    int w = image.cols, h = image.rows;
    int nr_superpixels = 400;
    int nc = 40;

    double step = sqrt((w * h) / (double) nr_superpixels);

    /* Perform the SLIC superpixel algorithm. */
    Slic slic;
    slic.generate_superpixels(lab_image, int(step), nc);
    slic.create_connectivity(lab_image);
    slic.display_contours(image, Vec3b(0, 0, 255));
    auto end_time = std::chrono::steady_clock::now();

    fout << "image: " << image_name << endl;
    fout << "size: " << w << "x" << h << " = " << w * h << " pixels" << endl;
    fout << "elapsed time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(
             end_time - start_time
         ).count() / 1000. << "s" << endl;;
    fout << "---------------------------------------------" << endl << endl;
    imwrite("out/" + image_name, image);

    /* Display the contours and show the result. */
//  imshow("result", image);
//  waitKey(0);
  }
  return 0;
}
