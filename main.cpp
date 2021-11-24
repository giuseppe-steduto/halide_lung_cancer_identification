// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"
//Include library to open/save images
#include "dcmtk/dcmimage/dicoimg.h"
using namespace Halide;

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

int main(int argc, char **argv) {
//    Buffer<uint8_t> input = Halide::Tools::load_image("SE0/1-001.dcm");

    //Halide::Tools::save_image(input, "output/1-001.dcm");
    std::cout << "Everything's okay!";
    return 0;
}
