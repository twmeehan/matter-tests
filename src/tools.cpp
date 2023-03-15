#include "tools.hpp"


unsigned int countlines(std::string file_name)
{
    std::ifstream file(file_name);

    std::string line;
    unsigned int p = 0; // particle

    if ( file.is_open() ) {
        while ( std::getline(file, line) ) {
            p++;
        }
    }
    else {
        std::cout << "Unable to open '"<< file_name << "'" << std::endl;
    }
    return p;
}

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
// linspace(a,b,N) return std::vector of length N in the closed interval [a,b], i.e., including b
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

// Taken from: https://stackoverflow.com/questions/21216909/these-python-functions-in-c
// Works like numpy.arange, does NOT include stop value
std::vector<T> arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}


bool copy_file(std::string source, std::string destination){
    std::ifstream in(source, std::ios::binary);
    std::ofstream out(destination, std::ios::binary);
    out << in.rdbuf();
    return in && out;
}
