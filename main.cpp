// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"
#include "tools/halide_image_io.h"
//Include library to open/save images
#include <dcmtk/dcmimgle/dcmimage.h>
#include <dcmtk/dcmdata/dcdatset.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcitem.h>

using namespace Halide;
using std::cout;
using std::cerr;
using std::endl;

Buffer<uint8_t> read_dicom_image(const char* filename) {
    DicomImage *image = new DicomImage(filename);
    if (image->getStatus() == EIS_Normal) {
        if (image->isMonochrome()) {
            image->setMinMaxWindow();
            auto *pixelData = (Uint8 *) (image->getOutputData(8 /* bits */));
            if (pixelData != nullptr) {
                Buffer<uint8_t> inp((uint8_t *) pixelData, 3, image->getWidth(), image->getHeight());
                return inp;
            }
            throw std::runtime_error("There is no data in the image!");
        }
    } else {
        cerr << "Error: cannot load DICOM image (" << DicomImage::getString(image->getStatus()) << ")"
                  << endl;
        delete image;
        throw std::runtime_error("Cannot load DICOM image!");
    }
}

Expr sobel (Func in)
{
    Func h("sobel_horiz"), v("sobel_vert");
    Func sob("sobel");
    Var x("x"), y("y"), c("c");
    h(x, y, c) = in(x+1, y-1, c) + 2*in(x+1, y, c) +
             in(x+1, y+1, c) - in(x-1, y-1, c) -
             2 * in(x-1, y, c) - in(x-1,y+1, c);
    v(x, y, c) = -in(x-1,y+1, c) - 2*in(x,y+1, c) -
             in(x+1,y+1, c) + in(x-1,y-1, c) +
             2 * in( x,y-1, c) + in( x+1,y-1, c);
    sob(x, y, c) = (abs(h(x,y,c)) +
                abs(v(x,y,c))) / 4;
    return sob(x, y, c);
}

Expr binarize (Func in) {
    Func bin("binarize");
    Var x("x"), y("y"), c("c");
    int threshold = 128;
    Expr tmp = in(x, y, c) - threshold;
    if (is_positive_const(tmp))
        return cast<uint8_t>(1);
    else
        return cast<uint8_t>(0);
}

int main(int argc, char **argv) {
    Buffer<uint8_t> image, complete;
    float gamma_exponent = 1.5; //Gamma correction constant
    Func gamma("gamma"), sobel_ed("sobel_edge_detecteor"), binarized("binarized_image");
    Var x("x"), y("y"), c("c"), xo("xo"), xi("xi"), yo("yo"), yi("yi");
    try {
        image = read_dicom_image("/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/1-001.dcm");
    } catch (const std::exception &e) {
        cerr << e.what() << endl;
    }
    gamma(x, y, c) = print(cast<uint8_t>(255 * pow(image(x, y, c) * 1.0f / 255, gamma_exponent)));
    sobel_ed(x, y, c) = sobel(gamma);
    binarized(x, y, c) = binarize(sobel_ed);

    complete = binarized.realize({image.width(), image.height(), image.channels()});

    //TODO other stuff

    //Create output file
    complete.copy_to_host();


    DcmFileFormat output_file;
    DcmDataset *dataset = output_file.getDataset();
    //TODO insert metadata
    DcmTag dcm_tag_pixeldata = DcmTag();
    Tools::convert_and_save_image(complete, "test.jpg");
    /*Uint8 image_buffer = save<uint8_t>(complete);
    dataset->putAndInsertUint8Array(dcm_tag_pixeldata, &complete, complete.size_in_bytes());
    OFCondition status = output_file.saveFile("test.dcm", EXS_LittleEndianExplicit);
    if (status.bad())
        cerr << "Error: cannot write DICOM file (" << status.text() << ")" << endl;*/

    return 0;
}
