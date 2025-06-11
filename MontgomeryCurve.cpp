#include "MontgomeryCurve.h"
#include <gmpxx.h>

MontgomeryCurve::MontgomeryCurve(const mpz_class& A, const mpz_class& n) : A(A), n(n) {
    mpz_class four = 4;
    mpz_class inv_four;
    if (mpz_invert(inv_four.get_mpz_t(), four.get_mpz_t(), n.get_mpz_t()) == 0) {
        throw std::runtime_error("Inverse of 4 mod n does not exist");
    }
    A24 = ((A + 2) * inv_four) % n;
}

MontgomeryPoint MontgomeryCurve::double_point(const MontgomeryPoint& P) const {
    mpz_class XpZ = (P.X + P.Z) % n;
    mpz_class XmZ = (P.X - P.Z + n) % n;

    mpz_class XpZ_sq; mpz_powm_ui(XpZ_sq.get_mpz_t(), XpZ.get_mpz_t(), 2, n.get_mpz_t());
    mpz_class XmZ_sq; mpz_powm_ui(XmZ_sq.get_mpz_t(), XmZ.get_mpz_t(), 2, n.get_mpz_t());

    mpz_class diff = (XpZ_sq - XmZ_sq + n) % n;

    mpz_class X2 = (XpZ_sq * XmZ_sq) % n;
    mpz_class Z2 = ((A24 * diff + XmZ_sq) % n * diff) % n;

    return MontgomeryPoint(X2, Z2);
}

MontgomeryPoint MontgomeryCurve::add_points(const MontgomeryPoint& P, const MontgomeryPoint& Q, const MontgomeryPoint& PminusQ) const {
    mpz_class XpZ = (P.X + P.Z) % n;
    mpz_class XmZ = (P.X - P.Z + n) % n;
    mpz_class XqZ = (Q.X + Q.Z) % n;
    mpz_class XqZm = (Q.X - Q.Z + n) % n;

    mpz_class t1 = (XpZ * XqZm) % n;
    mpz_class t2 = (XmZ * XqZ) % n;

    mpz_class X3 = ((t1 + t2) % n);
    mpz_class Z3 = ((t1 - t2 + n) % n);

    mpz_powm_ui(X3.get_mpz_t(), X3.get_mpz_t(), 2, n.get_mpz_t());
    mpz_powm_ui(Z3.get_mpz_t(), Z3.get_mpz_t(), 2, n.get_mpz_t());

    X3 = (X3 * PminusQ.X) % n;
    Z3 = (Z3 * PminusQ.Z) % n;

    return MontgomeryPoint(X3, Z3);
}

MontgomeryPoint MontgomeryCurve::scalar_multiply(const mpz_class& k, const MontgomeryPoint& P) const {
    MontgomeryPoint R0(1, 0);   // Point at infinity
    MontgomeryPoint R1 = P;

    size_t num_bits = mpz_sizeinbase(k.get_mpz_t(), 2);
    for (ssize_t i = static_cast<long>(num_bits) - 1; i >= 0; --i) {
        if (mpz_tstbit(k.get_mpz_t(), i) == 0) {
            R1 = add_points(R0, R1, P);
            R0 = double_point(R0);
        } else {
            R0 = add_points(R0, R1, P);
            R1 = double_point(R1);
        }
    }
    return R0;
}
