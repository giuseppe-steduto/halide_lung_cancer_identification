//
// Created by giuseppe on 18/12/21.
//

#ifndef LUNGCANCERIDENTIFICATION_ROI_DEFINITION_H
#define LUNGCANCERIDENTIFICATION_ROI_DEFINITION_H

// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"
#include "tools/halide_image_io.h"
//Include library to open/save images
#include <dcmtk/dcmimgle/dcmimage.h>
#include <dcmtk/dcmdata/dcdatset.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcitem.h>
#include "dcmtk/dcmdata/dctk.h"

using namespace Halide;
using std::cout;
using std::cerr;
using std::endl;

class DisjointSet
{
    std::unordered_map<int, int> parent;

    public:
        // perform MakeSet operation
        void makeSet(std::vector<int> const &universe);

        void add_element(int element);

        void add_equivalence(int a, int b);

        // find the root of the set in which element `k` belongs
        int find(int k);

        // Perform make_union of two subsets
        void make_union(int a, int b);

        bool is_root(int k);
};

int get_label(int a, int b, DisjointSet* ds);

Buffer<uint8_t> read_dicom_image(const char* filename);

Buffer<uint8_t> read_png_image(const char* filename);

Func mask (Func in, Buffer<uint8_t> mask);

Expr binarize (Func in, int threshold);

int get_otsu_treshold(Buffer<uint8_t> buffer);

void get_largest_cc(Buffer<uint8_t> buffer);

//TODO adjust background removal alogorithm
void background_removal(Buffer<uint8_t> mask);

int roi_definition(char* filename);

#endif //LUNGCANCERIDENTIFICATION_ROI_DEFINITION_H
