from misc import read_config, read_binary
import numpy as np
import cmath
import math
import os

NROW, NCOL = None, None
M_PI = math.pi
eps = np.finfo(np.float64).eps  # "machine epsilon" for 64-bit floating point number


t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p = None, None, None, None, None, None, None, None, None

def read_T3(d):
    global t11_p, t22_p, t33_p, t12_r_p, t12_i_p, t13_r_p, t13_i_p, t23_r_p, t23_i_p 
    sep = os.path.sep
    t11_p = read_binary(d + sep + 'T11.bin')
    t22_p = read_binary(d + sep + 'T22.bin')
    t33_p = read_binary(d + sep + 'T33.bin')
    t12_r_p = read_binary(d + sep + 'T12_real.bin')
    t12_i_p = read_binary(d + sep + 'T12_imag.bin')
    t13_r_p = read_binary(d + sep + 'T13_real.bin')
    t13_i_p = read_binary(d + sep + 'T13_imag.bin')
    t23_r_p = read_binary(d + sep + 'T23_real.bin')
    t23_i_p = read_binary(d + sep + 'T23_imag.bin')

def lamcloude(a, b, c, z1, z2, z3):
    p = 1./3.  # cf tra, z1p, z2p, z3p, fac0, fac1, fac2, fac3, s1, s2, deta, tr3;
    tra = (a + b + c) / 3.
    z1p, z2p, z3p = z1.conjugate(), z2.conjugate(), z3.conjugate() # conjugate(z1), conjugate(z2), conjugate(z3)

    fac0 = z1 * z1p + z2 * z2p + z3 * z3p

    s1 = a * b + a * c + b * c - fac0
    deta = a * b * c - c * z1 * z1p - b * z2 * z2p + z1 * z2p * z3 + z1p * z2 * z3p - a * z3 * z3p
    s2 = a * a - a * b + b * b - a * c - b * c + c * c + 3. * fac0
    fac1 = 27. * deta - 27. * s1 * tra + 54. * pow(tra, 3.)
    tr3 = fac1 + cmath.sqrt(pow(fac1, 2.) - 4. * pow(s2, 3.))
    fac2 = 1. + (1j * math.sqrt(3.))
    fac3 = 1. - (1j * math.sqrt(3.))

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
    n = math.sqrt(av1 * av1 + av2 * av2 + av3 * av3) + eps

    # normalised components as output
    v1, v2, v3 = v1 / n, v2 / n, v3 / n
    return [e1, e2, e3, v1, v2, v3]  #  double & e1, double & e2, double & e3, cf & v1, cf & v2, cf & v3


def rank1_t3(e1, v1, v2, v3): # e e1, cf v1, cf v2, cf v3):   #  generate T3 rank 1
    t11c = e1 * v1 * (v1.conjugate()) # conjugate(v1)
    t12c = e1 * v1 * (v2.conjugate()) # conjugate(v2)
    t13c = e1 * v1 * (v3.conjugate()) # conjugate(v3)
    t22c = e1 * v2 * (v2.conjugate()) # conjugate(v2)
    t23c = e1 * v2 * (v3.conjugate()) # conjugate(v3)
    t33c = e1 * v3 * (v3.conjugate()) # conjugate(v3)

    return [t11c, t12c, t13c, t22c, t23c, t33c]


def decom(i):   # calculate decom for pixel at linear index "i"
    # // intermediary variables
    # double t11, t12_r, t12_i, t13_r, t13_i, t22, t23_r, t23_i, t33;
    # double e1, e2, e3, p;
    # cf a, b, c, z1, z2, z3;
    t11 = t11_p[i]
    t22 = t22_p[i]
    t33 = t33_p[i]
    t12_r = t12_r_p[i]
    t12_i = t12_i_p[i]
    t13_r = t13_r_p[i]
    t13_i = t13_i_p[i]
    t23_r = t23_r_p[i]
    t23_i = t23_i_p[i]
    
    a = t11; b = t22; c = t33;
    z1 = t12_r + t12_i * 1j;
    z2 = t13_r + t13_i * 1j;
    z3 = t23_r + t23_i * 1j;
    
    # /* avoid 0 elements.. conditioning */
    eps2 = (a + b + c) * (1.0e-9) + eps;
    F = ((float(NROW) + float(NCOL)) / 2.);
    z1 = z1 + eps2 * F; # // %randn(sx,sy);
    z2 = z2 + eps2 * F; # // %randn(sx,sy);
    z3 = z3 + eps2 * F; #// %randn(sx,sy);
    a = a + eps2 * F; #// %randn(sx,sy);
    b = b + eps2 * F; #// %randn(sx,sy);
    c = c + eps2 * F; #// %randn(sx,sy);
    
    # //run lamcloude
    # cf v1, v2, v3;
    [e1, e2, e3, v1, v2, v3] = lamcloude(a, b, c, z1, z2, z3) #  e1, e2, e3, v1, v2, v3);
    
    # // rank 1 t3
    # cf t11c, t12c, t13c, t22c, t23c, t33c;
    [e1, e2, e3, v1, v2, v3] = rank1_t3(e1, v1, v2, v3) #  a, z1, z2, b, z3, c) #t11c, t12c, t13c, t22c, t23c, t33c);
    
    # // generate alpha etc. eigenvector parameters
    alpha = math.acos(abs(v1));
    phi = cmath.phase(t12c);
    theta = cmath.phase((t22c - t33c) + 2. * j * t23c.real) / 4.
    
    # // generate RGB colour composite from multiple eigenvector angles
    dn = alpha * 2. / M_PI # // alpha angle in red channel
    theta2 = theta + (theta > M_PI / 4.) * (M_PI / 4. - theta)
    theta2 = theta2 + (theta2 < -M_PI / 4.) * (-M_PI / 4. - theta2)
    vn = (theta2 + M_PI / 4.) * 2. / M_PI   # // az slope is green
    sn = abs(phi) / M_PI  # // mag of Pauli phase is blue (180 is Bragg)
    
    out_r[i] = dn;
    out_g[i] = vn;
    out_b[i] = sn;
    
    out_e1[i] = e1;
    out_e2[i] = e2;
    out_e3[i] = e3;
    
    #  // project data onto null channels // null_vecs=[o2d o3d];
    z1 = conjugate(o2d1)*v1 + conjugate(o2d2)*v2 + conjugate(o2d3)*v3;  # // oconj=o2d';
    z2 = conjugate(o3d1)*v1 + conjugate(o3d2)*v2 + conjugate(o3d3)*v3;  #// oconj=o3d';
    
    #  find optimum weights
    popt = cmath.phase(z2 * conjugate(z1)) * 180. / M_PI;
    za = (z1*conjugate(z1) - z2*conjugate(z2)) + j * 2.*abs(z1)*abs(z2);
    aopt = cmath.phase(za) * 90. / M_PI;
    ar = aopt * M_PI / 180.;
    br = popt * M_PI / 180.;
    
    # // optimum weight vector
    w1 = math.cos(ar) * o2d1 + math.sin(ar) * exp(j * br) * o3d1;
    w1 = conjugate(w1);
    w2 = math.cos(ar) * o2d2 + math.sin(ar) * exp(j * br) * o3d2;
    w2 = conjugate(w2);
    w3 = math.cos(ar) * o2d3 + math.sin(ar) * exp(j * br) * o3d3;
    w3 = conjugate(w3);
    
    # // find optimum subspace signal
    zopt = w1 * v1 + w2 * v2 + w3 * v3;
    ip = abs(zopt * conjugate(zopt));
    ip_eps = ip + eps;
    sopt = 10. * math.log(ip_eps) / math.log(10.); #// optimum normalised power
    
    sp = t11c + t22c + t33c; # //span power
    abs_sp = abs(sp);
    pwr = 10. * math.log(abs(sp)) / math.log(10.); # //span channel
    
    sm = math.fabs(t33);
    hv = 10. * math.log(sm) / math.log(10.); # log10(sm);
    sm2 = sopt + pwr;
    
    opt = pow(10., sm2 / 10.); #// linear opt channel
    
    out_opt[i] =  opt;
    out_v1[i] = sopt;
    '''
    //out_hv[i] = hv;
    //out_sm[i] = sm;
    
    //out_pwr[i] = (float)pwr;
    //out_sopt[i] = (float)sopt;
    //out_abs_sp[i] = (float)abs_sp;
    
    //out_ar[i] = (float) abs(ar);
    //out_br[i] = (float) abs(br);
    '''

x = read_config('../T3/config.txt')

NROW = x['nrow']
NCOL = x['ncol']

read_T3('../T3/')

decom(0)
