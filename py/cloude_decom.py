import numpy as np
import cmath

eps = np.finfo(np.float64).eps  # "machine epsilon" for 64-bit floating point number


def lamcloude(cf a, cf b, cf c, cf z1, cf z2, cf z3):
    p = 1./3.  # cf tra, z1p, z2p, z3p, fac0, fac1, fac2, fac3, s1, s2, deta, tr3;
    tra = (a + b + c) / 3.
    z1p, z2p, z3p = conjugate(z1), conjugate(z2), conjugate(z3)

    fac0 = z1 * z1p + z2 * z2p + z3 * z3p

    s1 = a * b + a * c + b * c - fac0
    deta = a * b * c - c * z1 * z1p - b * z2 * z2p + z1 * z2p * z3 + z1p * z2 * z3p - a * z3 * z3p
    s2 = a * a - a * b + b * b - a * c - b * c + c * c + 3. * fac0
    fac1 = 27. * deta - 27. * s1 * tra + 54. * pow(tra, 3.)
    tr3 = fac1 + cmath.sqrt(pow(fac1, 2.) - 4. * pow(s2, 3.))
    fac2 = 1. + j * sqrt(3.)
    fac3 = 1. - j * sqrt(3.)

    e1 = (tra + pow(tr3, p) / (3. * pow(2., p)) + (s2 * pow(2., p) + eps) / (3. * pow(tr3, p) + eps)).real
    e2 = (tra - (fac2 * s2) / (3. * pow(tr3, p) * pow(2., 2. * p) + eps) - (fac3 * pow(tr3,p)) / (6.*pow(2., p) + eps)).real
    e3 = (tra - (fac3 * s2) / (3. * pow(2., 2.*p) * pow(tr3, p) + eps) - (fac2 * pow(tr3,p)) / (6.*pow(2., p) + eps)).real
    
    if(e1 < e3):  # sort the eigenvalues
        tmp = e1; e1 = e3; e3 = tmp   #  swp(&e1, &e3); // sort eigenvalues
    if(e1 < e2):
        tmp = e1; e1 = e2; d2 = tmp   # swp(&e1, &e2);
    if(e2 < e3):
        tmp = e2; e2 = e3; e3 = tmp   # swp(&e2, &e3);
    if(e1 < e3 or e1 < e2 or e2 < e3):
        print("Warning: not sorted (%e, %e, %e)\n", e1, e2, e3);

    v2 = ((a - e1) * z3 - z1p * z2) / ((b - e1) * z2 - z3 * z1 + eps)  # dominant eigenvector
    v3 = (e1 - a - z1 * v2) / (z2 + eps)
    v1 = 1.  # ones(size(v2))

    av1, av2, av3 = abs(v1), abs(v2), abs(v3)
    n = sqrt(av1 * av1 + av2 * av2 + av3 * av3) + eps

    // normalised components as output
    v1, v2, v3 = v1 / n, v2 / n, v3 / n
    return [e1, e2, e3, v1, v2, v3]  #  double & e1, double & e2, double & e3, cf & v1, cf & v2, cf & v3


def rank1_t3(double e1, cf v1, cf v2, cf v3):   #  generate T3 rank 1
    t11c = e1 * v1 * conjugate(v1)
    t12c = e1 * v1 * conjugate(v2)
    t13c = e1 * v1 * conjugate(v3)
    t22c = e1 * v2 * conjugate(v2)
    t23c = e1 * v2 * conjugate(v3)
    t33c = e1 * v3 * conjugate(v3)

    return [t11c, t12c, t13c, t22c, t23c, t33c]
