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

Buffer<uint8_t> read_dicom_image(const char* filename) {
    DicomImage *image = new DicomImage(filename);
    if (image->getStatus() == EIS_Normal) {
        if (image->isMonochrome()) {
            image->setMinMaxWindow();
            uint8_t* pixelData = (uint8_t *) (image->getOutputData(8));
            if (pixelData != nullptr) {
                Buffer<uint8_t> inp((uint8_t *) pixelData, 3, image->getWidth(), image->getHeight());
                cout << "Width: " << image->getWidth() << " | Height: " << image->getHeight() << " | Size (bytes): " << inp.size_in_bytes() << endl;
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

Buffer<uint8_t> read_png_image(const char* filename) {
    return Tools::load_image(filename);
}

Func sobel (Func in)
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

Expr binarize (Func in, int threshold) {
    Func bin("binarize");
    Var x("x"), y("y");
    Expr tmp = in(x, y) - threshold;
    return select(in(x, y) > threshold, cast<uint8_t>(255), cast<uint8_t>(0));
}

int get_otsu_treshold(Buffer<uint8_t> buffer) {
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
    cout << "Thre: " << threshold << endl;
    return threshold;
}

int main(int argc, char **argv) {
    Buffer<uint8_t> image, complete;
    float gamma_exponent = 1.5; //Gamma correction constant
    Func gamma("gamma"), sobel_ed("sobel_edge_detecteor"), binarized("binarized_image"), eroded("eroded");
    Func h("sobel_horizontal"), v("sobel_vertical"), sobel_bounded("sobel_edge_detector_bounded");
    Var x("x"), y("y"), c("c"), xo("xo"), xi("xi"), yo("yo"), yi("yi");
    try {
//        image = read_dicom_image("/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/1-001.dcm");
        image = read_png_image("/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/1-001.jpeg");
    } catch (const std::exception &e) {
        cerr << "Errore nell'apertura del file:" << endl;
        cerr << e.what() << endl;
    }

    //Gamma correction
    gamma(x, y) = cast<uint8_t>(255 * pow(image(x, y) * 1.0f / 255, gamma_exponent));

    //Sobel edge detection
//    //TODO border handling for sobel edge detection
//    h(x, y) = gamma(x + 1, y - 1) + 2 * gamma(x + 1, y) +
//                 gamma(x + 1, y + 1) - gamma(x - 1, y - 1) -
//                 2 * gamma(x - 1, y) - gamma(x - 1, y + 1);
//    v(x, y) = -gamma(x - 1, y + 1) - 2 * gamma(x, y + 1) -
//                 gamma(x + 1, y + 1) + gamma(x - 1, y - 1) +
//                 2 * gamma(x, y - 1) + gamma(x + 1, y - 1);
//    sobel_ed(x, y) = (h(x, y) + v(x, y)) / 4;
    //Boundary conditions will also erode the image
    //sobel_bounded(x, y) = BoundaryConditions::constant_exterior(sobel_ed, erased_value, 2, image.width() - 4, 2, image.height() - 4)(x, y);

    //sobel_bounded.trace_stores();
    //Buffer<uint8_t> tmp(image.width(), image.height(), image.channels());
    //sobel_bounded.realize(tmp);

    //Image binarization and erosion
    int threshold = get_otsu_treshold(image);
    binarized(x, y) = cast<uint8_t>(binarize(gamma, threshold));
    uint8_t erased_value = 0;
    eroded(x, y) = BoundaryConditions::constant_exterior(binarized, erased_value, 2, image.width() - 4, 2, image.height() - 4)(x, y);

    //TODO get largest connected component

    //TODO complement image (to be used as mask)

    //TODO mask image

    complete = eroded.realize({image.width(), image.height(), image.channels()});

    Tools::save_image(complete, "/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/dcm_modified.jpeg");
    return 0;

    //Output file
    DcmFileFormat output_file;
    DcmDataset *dataset = output_file.getDataset();
    char uid[100];
    //TODO insert metadata

    dataset->putAndInsertUint8Array(DCM_PixelData, image.begin(), image.size_in_bytes());
    dataset->putAndInsertString(DCM_SOPClassUID, UID_SecondaryCaptureImageStorage);
    dataset->putAndInsertString(DCM_SOPInstanceUID, dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT));
//    dataset->putAndInsertUint8Array(DCM_PixelData, complete.begin(), complete.size_in_bytes());
    OFCondition status = output_file.saveFile("/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/test.dcm", EXS_BigEndianExplicit);
    if (status.bad())
        cerr << "Error: cannot write DICOM file (" << status.text() << ")" << endl;
    else
        cout << "I should have written the DICOM file :)" << endl;

    return 0;
}
