#include "tools.hpp"

unsigned int load_array(std::vector<TV>& array, std::string file_name)
{
    std::ifstream file(file_name);

    std::string line;
    T value;
    unsigned int p = 0; // particle
    unsigned int j;     // component (x, y or z)

    if ( file.is_open() ) {
        while ( std::getline(file, line) ) {
            j = 0;
            std::stringstream line_stream(line);
            while ( line_stream >> value ) {
                // debug("p = ", p, " j = ", j, " value = ", value);
                array[p](j) = value;
                j++;
            }
            p++;
        }
    }
    else {
        std::cout << "Unable to open '"<< file_name << "'" << std::endl;
    }

    return p;
}


// Taken from: https://gist.github.com/lorenzoriano/5414671
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}
