#include <filesystem>
#include <chrono>
#include "roi_definition.cpp"

#define PATH "/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/SE1/"
#define PATH_OUTPUT "/home/giuseppe/Desktop/PoliMi/NECSTCamp/LungCancerIdentification/output/"

int main(int argc, char **argv) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (const auto & entry : std::filesystem::directory_iterator(PATH)) {
        char filename[255];
        strcpy(filename, entry.path().string().c_str());
        roi_definition((char *) filename, PATH_OUTPUT);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
}
