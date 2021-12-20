// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"
#include "roi_definition.h"
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

inline void DisjointSet::makeSet(std::vector<int> const &universe)
{
    // create `n` disjoint sets (one for each item)
    for (int i: universe) {
        parent[i] = i;
    }
}

inline void DisjointSet::add_element(int element) {
    parent[element] = element;
}

inline void DisjointSet::add_equivalence(int a, int b) {
    try {
        int tmp = parent.at(a);
    } catch (const std::out_of_range& e) {
        add_element(a);
    }
    try {
        int tmp = parent.at(b);
    } catch (const std::out_of_range& e) {
        add_element(b);
    }
    make_union(a, b);
}

// find the root of the set in which element `k` belongs
inline int DisjointSet::find(int k)
{
    // if `k` is root
    if (parent[k] == k) {
        return k;
    }

    // recur for the parent until we find the root
    return find(parent[k]);
}

// Perform make_union of two subsets
inline void DisjointSet::make_union(int a, int b)
{
    // find the root of the sets in which elements `x` and `y` belongs
    if (a == 0 || b == 0)
        return;
    int x = find(a);
    int y = find(b);

    parent[x] = y;
}

inline bool DisjointSet::is_root(int k) {
    if (parent[k] == k)
        return true;
    return false;
}

inline Buffer<uint8_t> read_dicom_image(const char* filename) {
    DicomImage *image = new DicomImage(filename);
    if (image->getStatus() == EIS_Normal) {
        if (image->isMonochrome()) {
            image->setMinMaxWindow();
            uint8_t* pixelData = (uint8_t *) (image->getOutputData(8));
            if (pixelData != nullptr) {
                Buffer<uint8_t> inp((uint8_t *) pixelData, image->getWidth(), image->getHeight(), 1);
                cout << "Width: " << image->getWidth() << " | Height: " << image->getHeight() << " | Size (bytes): " << inp.size_in_bytes() << " --- ";
                return inp;
            }
            throw std::runtime_error("There is no data in the image!");
        } else {
            cerr << "Image is not monochrome!" << endl;
        }
    } else {
        cerr << "Error: cannot load DICOM image (" << DicomImage::getString(image->getStatus()) << ")"
             << endl;
        delete image;
        throw std::runtime_error("Cannot load DICOM image!");
    }
}

inline Buffer<uint8_t> read_png_image(const char* filename) {
    return Tools::load_image(filename);
}

inline Func sobel (Func in)
{
    Func h("sobel_horiz"), v("sobel_vert");
    Func sob("sobel");
    Var x("x"), y("y"), c("c");

    h(x, y, c) = in(x + 1, y - 1, c) + 2 * in(x + 1, y, c) +
                 in(x + 1, y + 1, c) - in(x - 1, y - 1, c) -
                 2 * in(x - 1, y, c) - in(x - 1, y + 1, c);
    v(x, y, c) = -in(x - 1, y + 1, c) - 2 * in(x, y + 1, c) -
                 in(x + 1, y + 1, c) + in(x - 1, y - 1, c) +
                 2 * in(x, y - 1, c) + in(x + 1, y - 1, c);
    sob(x, y, c) = (h(x, y, c) +
                    v(x, y, c)) / 4;
    return sob;
}

inline Func mask (Func in, Buffer<uint8_t> mask) {
    Var x("x"), y("y"), c("c");
    Func mask_func("mask");
    mask_func(x, y, c) = select(mask(x, y) == 255, in(x, y, c), 0);
    return mask_func;
}

inline Expr binarize (Func in, int threshold) {
    Func bin("binarize");
    Var x("x"), y("y"), c("c");
    Expr tmp = in(x, y, c) - threshold;
    return select(in(x, y, c) > threshold, cast<uint8_t>(255), cast<uint8_t>(0));
}

inline int get_otsu_treshold(Buffer<uint8_t> buffer) {
    int histogram[256] = {0};

    //Fill in histogram
    for (int x = 0; x < buffer.width(); x++) {
        for (int y = 0; y < buffer.height(); y++) {
            histogram[buffer(x, y)]++;
        }
    }

    int total = buffer.size_in_bytes(); //Total number of pixels

    float sum = 0;
    for (int t=0 ; t<256 ; t++) {
        sum += t * histogram[t];
    }

    float sumB = 0;
    int wB = 0;
    int wF = 0;

    float varMax = 0;
    int threshold = 0;

    for (int t=0 ; t<256 ; t++) {
        wB += histogram[t];               // Weight Background
        if (wB == 0) continue;
        wF = total - wB;                 // Weight Foreground
        if (wF == 0) continue;
        sumB += (float) (t * histogram[t]);
        float mB = sumB / wB;            // Mean Background
        float mF = (sum - sumB) / wF;    // Mean Foreground

        // Calculate Between Class Variance
        float varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);

        // Check if new maximum found
        if (varBetween >= varMax) {
            varMax = varBetween;
            threshold = t;
        }
    }
    cout << "Otsu threshold: " << threshold << " --- ";
    return threshold;
}

inline void get_largest_cc(Buffer<uint8_t> buffer) {
    DisjointSet equivalences;
    int map[buffer.width()][buffer.height()]; //Matrix of labels for the pixels
    std::unordered_map<int, int> labels; //How many pixels for each label
    std::pair<int, int> max = std::make_pair(0, 0);

    //First step: assign label to every pixel (0 for bg) and annote equivalences
    for (int i = 0; i < buffer.width(); i++)
        map[0][i] = 0;
    for (int j = 0; j < buffer.height(); j++)
        map[j][0] = 0;

    int last_label = 0;
    for (int x = 1; x < buffer.width(); x++) {
        for (int y = 1; y < buffer.height() - 1; y++) {
            //If pixel is background, skip to the next one
            if (buffer(x, y) == 0) {
                map[x][y] = 0;
                continue;
            }

            //Assign label to pixel buffer(x, y) by checking the one above and the one to the left
            if (map[x-1][y] == 0 && map[x][y-1] == 0) {
                //Assign new label to pixel
                last_label++;
                map[x][y] = last_label;
            }
            else {
                int selected_label = get_label(map[x-1][y], map[x][y-1], &equivalences);
                map[x][y] = selected_label;
            }
        }
    }

    //Now take care of the equivalences
    for (int x = 1; x < buffer.width(); x++) {
        for (int y = 1; y < buffer.height() - 1; y++) {
            if (map[x][y] == 0)
                continue;
            if (!equivalences.is_root(map[x][y])) {
                map[x][y] = equivalences.find(map[x][y]);
            }
            try {
                labels.at(map[x][y])++;
            } catch (const std::out_of_range& e) {
                labels[map[x][y]] = 1;
            }
            if (labels.at(map[x][y]) > max.second) {
                max.first = map[x][y];
                max.second = labels.at(map[x][y]);
            }
        }
    }

    //Now modify image to let only largest component remain
    for (int x = 1; x < buffer.width(); x++) {
        for (int y = 1; y < buffer.height() - 1; y++) {
            if (map[x][y] != max.first) {
                buffer(x, y) = 0;
            } else {
                buffer(x, y) = 255;
            }
        }
    }
}

//TODO adjust background removal alogorithm
inline void background_removal(Buffer<uint8_t> mask) {
    int map[mask.width()][mask.height()]; //Matrix of labels for the pixels

    bool isBg;
    for (int i = 0; i < mask.height(); i++) {
        int prev = mask(i, 0);
        //From the beginning of the row (border), till another 0 is met (thorax), make
        //all pixels 0 (as they're part of the background)
        isBg = true;
        for (int j = 1; isBg && j < mask.width(); j++) {
            int curr = mask(i, j);
            if (curr == 0) { //Is already removed from mask
                prev = curr;
                continue;
            }

            if (curr == prev) {
                prev = curr;
                mask(i, j) = 0;
            } else {
                isBg = false;
            }
        }

        //Make the same calculation from the end of the row
        isBg = true;
        prev = mask(i, mask.width() - 1);
        for(int j = mask.width() - 1; isBg && j >= 1; j--) {
            int curr = mask(i, j);
            if (curr == 0) { //Is already removed from mask
                prev = 0;
                continue;
            }

            if (curr == prev) {
                prev = curr;
                mask(i, j) = 0;
            } else {
                isBg = false;
            }
        }
    }
}

inline int get_label(int a, int b, DisjointSet* ds) {
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == b) return a;

    //a and b have different labels, different from 0
    ds->add_equivalence(a, b);
    return a;
}

inline int roi_definition(char* filename) {
    char* file_path;
    Buffer<uint8_t> image, mask_buff, masked_buff;
    float gamma_exponent = 1.5; //Gamma correction constant
    Func gamma("gamma"), sobel_ed("sobel_edge_detecteor"), binarized("binarized_image"), eroded("eroded");
    Func masked("masked"), complement("complement");
    Func h("sobel_horizontal"), v("sobel_vertical"), sobel_bounded("sobel_edge_detector_bounded");
    Var x("x"), y("y"), c("c"), xo("xo"), xi("xi"), yo("yo"), yi("yi");
    try {;
        file_path = filename;
        image = read_dicom_image(file_path);
    } catch (const std::exception &e) {
        cerr << "Errore nell'apertura del file:" << endl;
        cerr << e.what() << endl;
        return 1;
    }

    Tools::save_image(image, strcat(file_path, "converted.jpg"));

    //Gamma correction
    gamma(x, y, c) = cast<uint8_t>(255 * pow(image(x, y, c) * 1.0f / 255, gamma_exponent));

    //Image binarization
    int threshold = get_otsu_treshold(image);
    binarized(x, y, c) = cast<uint8_t>(binarize(gamma, threshold));
    uint8_t erased_value = 0;
    eroded(x, y, c) = BoundaryConditions::constant_exterior(binarized, erased_value, 2, image.width() - 4, 2, image.height() - 4)(x, y, c);

    //Realize mask buffer and mask image
    mask_buff = eroded.realize({image.width(), image.height(), image.channels()});
    get_largest_cc(mask_buff);

    //Complement mask
    complement(x, y, c) = cast<uint8_t>(select(mask_buff(x, y, c) == 255, 0, 255));

    //Background removal from complemented mask
    mask_buff = complement.realize({image.width(), image.height(), image.channels()});
    background_removal(mask_buff);
    Tools::save_image(mask_buff, strcat(file_path, "tmp_mask.jpg")); //TODO REMOVE THIS

    masked(x, y, c) = select(mask_buff(x, y, c) == 255, gamma(x, y, c), 0);
    masked_buff = masked.realize({image.width(), image.height(), image.channels()});

    //Save masked image
    Tools::save_image(masked_buff, strcat(file_path, "dcm_masked.jpg"));
//    return 0;

    //Output file
    DcmFileFormat output_file;
    DcmDataset *dataset = output_file.getDataset();
    char uid[100];
    //TODO insert metadata

    dataset->putAndInsertUint8Array(DCM_PixelData, image.begin(), image.size_in_bytes());
    dataset->putAndInsertString(DCM_SOPClassUID, UID_SecondaryCaptureImageStorage);
    dataset->putAndInsertString(DCM_SOPInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
//    dataset->putAndInsertUint8Array(DCM_PixelData, mask_buff.begin(), mask_buff.size_in_bytes());
    OFCondition status = output_file.saveFile(strcat(file_path, "test.dcm"), EXS_BigEndianExplicit);
    if (status.bad())
        cerr << "Error: cannot write DICOM file (" << status.text() << ")" << endl;
    else
        cout << "I should have written the DICOM file :)" << endl;

    return 0;
}