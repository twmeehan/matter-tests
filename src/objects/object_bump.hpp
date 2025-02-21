// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTBUMP_HPP
#define OBJECTBUMP_HPP

#include "object_general.hpp"

class ObjectBump : public ObjectGeneral {
public:

    ~ObjectBump(){}

    ObjectBump(BoundaryCondition bc_in, T friction_in, std::string name_in = "") : ObjectGeneral(bc_in, friction_in, name_in) {}

    bool inside(const TV& X_in) override {

        T x = X_in(0);
        T y = X_in(1);

        T y_limit = 0.0475 / std::cosh(25*(x-0.43));

        if (y < y_limit)
            return true;
        else
            return false;

    }

    TV normal(const TV& X_in) override {

        T x = X_in(0);
        TV n;

        T arg = 25*(x-0.43);
        T b_der = -25 * 0.0475 * std::tanh(arg) / std::cosh(arg);
        n(0) = -b_der;
        n(1) = 1;

        return n.normalized();
    }

};

#endif  // OBJECTBUMP_HPP
