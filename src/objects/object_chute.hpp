#ifndef OBJECTCHUTE_HPP
#define OBJECTCHUTE_HPP

#include "object_general.hpp"

class ObjectChute : public ObjectGeneral {
public:
    T l1;
    T l2;
    T l3;
    T a1;
    T a2;
    T a3;

    ~ObjectChute(){}

    ObjectChute(BoundaryCondition bc_in, T friction_in, std::string name_in, T l1_in = 20, T l2_in = 4.2, T l3_in = 6.05, T a1_in = 13, T a2_in = 20, T a3_in = 37) : ObjectGeneral(bc_in, friction_in, name_in), l1(l1_in), l2(l2_in), l3(l3_in), a1(a1_in*M_PI/180.0), a2(a2_in*M_PI/180.0), a3(a3_in*M_PI/180.0)  {}

    bool inside(const TV& X_in) override {

        T x = X_in(0);
        T y = X_in(1);

        T L1 = l1;
        T L2 = l2*std::cos(a1);
        T L3 = l3*std::cos(a2);

        if ( (x < L1)            && (y < 0) )
            return true;
        else if( (x < L1+L2)     && (y < std::tan(a1)*(x-L1)) )
            return true;
        else if( (x < L1+L2+L3)  && (y < std::tan(a2)*(x-L1-L2) + l2*std::sin(a1)) )
            return true;
        else if( (x >= L1+L2+L3) && (y < std::tan(a3)*(x-L1-L2-L3) + l2*std::sin(a1) + l3*std::sin(a2)) )
            return true;
        else
            return false;

    }

    TV normal(const TV& X_in) override {

        T x = X_in(0);
        T y = X_in(1);
        TV n;


        T L1 = l1;
        T L2 = l2*std::cos(a1);
        T L3 = l3*std::cos(a2);

        if ( (x < L1)            && (y < 0) ){
            n(0) = 0;
            n(1) = 1;
        }
        else if( (x < L1+L2)     && (y < std::tan(a1)*(x-L1)) ){
            n(0) = -std::tan(a1);
            n(1) = 1;
        }
        else if( (x < L1+L2+L3)  && (y < std::tan(a2)*(x-L1-L2) + l2*std::sin(a1)) ){
            n(0) = -std::tan(a2);
            n(1) = 1;
        }
        else if( (x >= L1+L2+L3) && (y < std::tan(a3)*(x-L1-L2-L3) + l2*std::sin(a1) + l3*std::sin(a2)) ){
            n(0) = -std::tan(a3);
            n(1) = 1;
        }
        else{ // should not occur
            n(0) = 0;
            n(1) = 0;
        }

        return n.normalized();
    }

};

#endif  // OBJECTCHUTE_HPP
