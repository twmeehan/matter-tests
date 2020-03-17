#include "tools.hpp"

void load_array(TVX& array_, unsigned int n_cols, std::string file_name)
{
    std::ifstream file_(file_name);

    if ( file_.is_open() ) {

        std::string line_;
        T value_;
        unsigned int i = 0, j;
        while ( std::getline(file_, line_) ) {

            j = 0;
            std::stringstream line_stream(line_);
            while ( line_stream >> value_ ) {
                array_[i*n_cols + j] = value_;
                j++;
            }
            i++;
        }
    }
    else {
        std::cout << "Unable to open '"<< file_name << "'" << std::endl;
    }
}
