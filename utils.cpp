#include<random>

std::mt19937_64 utils::get_RNG()
{
    std::random_device rd;
    std::mt19937_64 g(rd());
    return g;
} 
