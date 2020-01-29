#include <iostream>
#include "dynamics.hpp"

int main()
{
    pcg32 RNG(23, 42);
    struct trial trial;
    initTrial(trial, 10, 2, 2.0, RNG, NULL, NULL);
    std::cout << "MAXN=" << MAXN << std::endl;
    std::cout << "MAXK=" << MAXK << std::endl;
    std::cout << "sizeof(trial)[MB]=" << sizeof(trial)/1e3 << std::endl;
    return 0;
}
