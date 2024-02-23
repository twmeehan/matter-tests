#ifndef GENERALOBJECT_HPP
#define GENERALOBJECT_HPP

#include "tools.hpp"
#include <string>


class GeneralObj{
public:

    BoundaryCondition bc;
    T friction;
    std::string name;

    GeneralObj(BoundaryCondition bc, T friction, std::string name) : bc(bc), friction(friction), name(name) {}

    virtual ~GeneralObj(){}

    virtual bool inside(const TV& X_in) = 0;

    virtual TV normal(const TV& X_in) = 0;

};

#endif  // GENERALOBJECT_HPP
