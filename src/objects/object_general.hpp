// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTGENERAL_HPP
#define OBJECTGENERAL_HPP

#include "../tools.hpp"

class ObjectGeneral{
public:

    BoundaryCondition bc;
    T friction;
    std::string name;

    ObjectGeneral(BoundaryCondition bc, T friction, std::string name) : bc(bc), friction(friction), name(name) {}

    virtual ~ObjectGeneral(){}

    virtual bool inside(const TV& X_in) = 0;

    virtual TV normal(const TV& X_in) = 0;

};

#endif  // OBJECTGENERAL_HPP
