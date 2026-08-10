# -*- coding: utf-8 -*-
"""Eshelby tensor for a spheroid, by numerical quadrature, and the dilute
estimate for a composite of ALIGNED spheroids.

Round-1 benchmarked the isotropic knockdown against spherical Mori-Tanaka. The
paper's headline claim is an anisotropy produced by aligned, non-spherical
pockets, and spheres cannot test that. This computes the effective transversely
isotropic tensor for aligned spheroids and predicts E_z/E_x, which can be put
against the measured value.

Validated against the closed-form spherical Eshelby components before use.
"""
import numpy as np

# ---- matrix
E_m, nu_m = 9.4, 0.33
K_m = E_m / (3 * (1 - 2 * nu_m))
G_m = E_m / (2 * (1 + nu_m))
lam = K_m - 2.0 / 3.0 * G_m


def iso_C(lmbda, mu):
    C = np.zeros((3, 3, 3, 3))
    d = np.eye(3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i, j, k, l] = (lmbda * d[i, j] * d[k, l]
                                     + mu * (d[i, k] * d[j, l] + d[i, l] * d[j, k]))
    return C


Cm = iso_C(lam, G_m)


def eshelby(a, C, n_theta=120, n_phi=240):
    """Eshelby tensor for an ellipsoid with semi-axes a=(a1,a2,a3) in matrix C.

    Gavazzi & Lagoudas quadrature: integrate over the unit sphere, mapping each
    direction through the semi-axes.
    """
    S = np.zeros((3, 3, 3, 3))
    # Gauss-Legendre in cos(theta), uniform in phi
    x, wx = np.polynomial.legendre.leggauss(n_theta)
    phis = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    dphi = 2 * np.pi / n_phi
    a = np.asarray(a, float)
    for ct, w in zip(x, wx):
        st = np.sqrt(max(0.0, 1 - ct * ct))
        for ph in phis:
            n = np.array([st * np.cos(ph), st * np.sin(ph), ct])
            xi = n / a
            # acoustic tensor K_ik = C_ijkl xi_j xi_l
            K = np.einsum('ijkl,j,l->ik', C, xi, xi)
            D = np.linalg.inv(K)
            # G_ijkl = D_ik xi_j xi_l  (symmetrised on ij and kl below)
            G = np.einsum('ik,j,l->ijkl', D, xi, xi)
            Gs = 0.25 * (G + G.transpose(1, 0, 2, 3)
                         + G.transpose(0, 1, 3, 2) + G.transpose(1, 0, 3, 2))
            S += w * dphi * np.einsum('ijmn,mnkl->ijkl', Gs, C)
    return S / (4 * np.pi)


def to_voigt(T):
    m = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    V = np.zeros((6, 6))
    for I, (i, j) in enumerate(m):
        for J, (k, l) in enumerate(m):
            V[I, J] = T[i, j, k, l] * (2 if J >= 3 else 1)
    return V


# ---------------------------------------------------------------- validation
S_sph = eshelby((1.0, 1.0, 1.0), Cm)
exact_1111 = (7 - 5 * nu_m) / (15 * (1 - nu_m))
exact_1122 = (5 * nu_m - 1) / (15 * (1 - nu_m))
print('validation against the closed-form spherical Eshelby tensor (nu=%.2f)' % nu_m)
print('   S_1111  numeric %.5f   exact %.5f   err %.1e'
      % (S_sph[0, 0, 0, 0], exact_1111, abs(S_sph[0, 0, 0, 0] - exact_1111)))
print('   S_1122  numeric %.5f   exact %.5f   err %.1e'
      % (S_sph[0, 0, 1, 1], exact_1122, abs(S_sph[0, 0, 1, 1] - exact_1122)))
print()


def dilute_aligned(aspect, K_i, G_i, phi):
    """Dilute estimate for aligned spheroids, semi-axes (1,1,aspect) along z."""
    Ci = iso_C(K_i - 2.0 / 3.0 * G_i, G_i)
    S = eshelby((1.0, 1.0, aspect), Cm)
    Sv, Cmv, Civ = to_voigt(S), to_voigt(Cm), to_voigt(Ci)
    dC = Civ - Cmv
    # strain concentration A = [I + S Cm^-1 (Ci - Cm)]^-1
    A = np.linalg.inv(np.eye(6) + Sv @ np.linalg.inv(Cmv) @ dC)
    return Cmv + phi * dC @ A


def moduli(Cv):
    Sc = np.linalg.inv(Cv)
    return 1.0 / Sc[0, 0], 1.0 / Sc[2, 2]      # E_x, E_z


BRINE = (2.2, 0.44e-3)
print('aligned spheroids along z: predicted anisotropy at the column brine fraction')
print('%10s %10s %10s %10s' % ('aspect', 'E_x GPa', 'E_z GPa', 'E_z/E_x'))
for asp in (1.0, 1.3, 1.6, 2.0, 2.5):
    Cv = dilute_aligned(asp, *BRINE, 0.040)
    Ex, Ez = moduli(Cv)
    print('%10.2f %10.3f %10.3f %10.4f' % (asp, Ex, Ez, Ez / Ex))
print()
print('and horizontally aligned (spheroid axis along x), same fraction:')
for asp in (1.6, 2.0):
    S = eshelby((asp, 1.0, 1.0), Cm)
    Ci = iso_C(BRINE[0] - 2.0 / 3.0 * BRINE[1], BRINE[1])
    Sv, Cmv, Civ = to_voigt(S), to_voigt(Cm), to_voigt(Ci)
    dC = Civ - Cmv
    A = np.linalg.inv(np.eye(6) + Sv @ np.linalg.inv(Cmv) @ dC)
    Cv = Cmv + 0.040 * dC @ A
    Ex, Ez = moduli(Cv)
    print('%10.2f %10.3f %10.3f %10.4f' % (asp, Ex, Ez, Ez / Ex))
