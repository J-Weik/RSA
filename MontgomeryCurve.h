//
// Created by Justus Weik on 11.06.2025.
//

#ifndef MONTGOMERYCURVE_H
#define MONTGOMERYCURVE_H
#include <gmpxx.h>

struct MontgomeryPoint {
    mpz_class X;
    mpz_class Z;
};
class MontgomeryCurve {
public:
    MontgomeryCurve(const mpz_class& A, const mpz_class& n);
    [[nodiscard]] MontgomeryPoint double_point(const MontgomeryPoint& P) const;
    [[nodiscard]] MontgomeryPoint add_points(const MontgomeryPoint& P, const MontgomeryPoint& Q, const MontgomeryPoint& PminusQ) const;
    [[nodiscard]] MontgomeryPoint scalar_multiply(const mpz_class& k, const MontgomeryPoint& P) const;
private:
    mpz_class A;  // Kurvenparameter A
    mpz_class n;  // Modul n
    mpz_class A24;
};

#endif //MONTGOMERYCURVE_H
