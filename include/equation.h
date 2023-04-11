#ifndef LAB2_V2_INCLUDE_EQUATION_H_
#define LAB2_V2_INCLUDE_EQUATION_H_

template<unsigned Order>
struct LinearODE {;
    double coeffs[Order + 2];

    inline unsigned order() const { return Order; }
};



#endif //LAB2_V2_INCLUDE_EQUATION_H_
