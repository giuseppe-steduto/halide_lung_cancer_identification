// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"
//Include library to open/save images
#include <dcmtk/dcmimgle/dcmimage.h>

using namespace Halide;
using std::cout;

//Gamma correction
void gamma(const uint8_t* src, uint8_t* dst, int height, int width, float gamma, bool gpu) {
    static Func f("gamma_correction");
    static Buffer<uint8_t> inp((uint8_t*) src, 3, width, height), out(dst, 3, width, height);

    if (!f.defined()) {
        Var x("x"), y("y"), c("c"), xo("xo"), xi("xi"), yo("yo"), yi("yi");
        f(c, x, y) = cast<uint8_t>(255 * pow(inp(c, x, y) * 1.0f / 255, gamma));
        Target t = Halide::get_host_target();
        if (gpu) {
            t.set_feature(Halide::Target::OpenCL);
            f.bound(x, 0, width).bound(y, 0, height).bound(c, 0, 3)
                    .split(x, xo, xi, 16).split(y, yo, yi, 16).reorder(xi, yi, c, xo, yo)
                    .gpu_blocks(c, xo, yo).gpu_threads(xi, yi);
        } else {
            f.set_estimate(x, 0, width).set_estimate(y, 0, height).set_estimate(c, 0, 3);
            Pipeline(f).auto_schedule(t);
        }
        f.compile_jit(t);
    }
    inp.set_host_dirty();
    f.realize(out);
    out.copy_to_host();
}

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
        std::cerr << "Error: cannot load DICOM image (" << DicomImage::getString(image->getStatus()) << ")"
                  << std::endl;
        delete image;
        throw std::runtime_error("Cannot load DICOM image!");
    }
}

int main(int argc, char **argv) {
    Buffer<uint8_t> image;
    float gamma_exponent = 1.5; //Gamma correction constant
    Func gamma;
    Var x("x"), y("y"), c("c"), xo("xo"), xi("xi"), yo("yo"), yi("yi");
    try {
        image = read_dicom_image("/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/1-001.dcm");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    gamma(c, x, y) = cast<uint8_t>(255 * pow(image(c, x, y) * 1.0f / 255, gamma_exponent));
    gamma.realize();



    return 0;
}
