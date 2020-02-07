#include <cstdio>
#include <string>
#include <time.h>
#include <fstream>

std::string defaultFilePath(int N, int K, const char* affix, const char* suffix)
{
    char filename[128];
    char path[128];
    char date[64];
    time_t tt;
    time(&tt);
    struct tm* ti = localtime(&tt);
    strftime(date, 64, "%F-%b-%a", ti);

    int version = 0;
    sprintf(path, "data/%s/", date);
    std::string cmd = "mkdir -p " + std::string(path);
    system(cmd.c_str());

    sprintf(filename, "%s-N-%06dK-%06d%s_v%d.dat", affix, N, K, suffix, version);
    while (std::ifstream(std::string(path) + std::string(filename))) {
        version++;
        sprintf(filename, "%s-N-%06dK-%06d%s_v%d.dat", affix, N, K, suffix, version);
    }

    return std::string(path) + std::string(filename);
}
