#include <filesystem>
#include <chrono>
#include "roi_definition.cpp"

#define PATH "example/input/path/folder" //Change this to the directory where your scans are located
#define PATH_OUTPUT "example/output/path/folder/" //Change this to the directory where you want to output the segmented images

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
