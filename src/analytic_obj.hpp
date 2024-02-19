#ifndef ANALYTICOBJECT_HPP
#define ANALYTICOBJECT_HPP

#include "tools.hpp"
#include <string>


class AnalyticObj{
public:

    AnalyticObj(BoundaryCondition bc, T friction, std::string name, int type, T h = 0.016) : bc(bc), friction(friction), name(name), type(type), h(h) {}

    bool inside(const TV& X_in){

        T x = X_in(0);
        T y = X_in(1);

        if (type == 0){     // viroulet bump
            T y_limit = 0.0475 / std::cosh(25*(x-h));

            if (y < y_limit)
                return true;
            else
                return false;

        } else if (type > 0){ // quad gate
            T y_limit = h + type * x*x;

            if (y > y_limit)
                return true;
            else
                return false;

        } else{ // type < 0   // ramp
            T y_limit = 0.1 * std::tanh(type*x);
            if (y < y_limit)
                return true;
            else
                return false;

        }
    }

    TV normal(const TV& X_in){

        T x = X_in(0);
        TV n;

        if (type == 0){
            T arg = 25*(x-h);
            T b_der = -25 * 0.0475 * std::tanh(arg) / std::cosh(arg);
            n(0) = -b_der;
            n(1) = 1;
        } else if (type > 0){
            T b_der = 2.0 * type * x;
            n(0) = -b_der;
            n(1) = 1;
            n *= -1;
        } else{ // type < 0
            T tmp = std::cosh(type*x);
            T b_der = 0.1*type/(tmp*tmp);
            n(0) = -b_der;
            n(1) = 1;
        }

        return n.normalized();
    }

    BoundaryCondition bc;
    T friction;
    std::string name;
    int type;
    T h;

};

#endif  // ANALYTICOBJECT_HPP
