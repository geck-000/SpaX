# -*- coding: utf-8 -*-
"""Benchmark the measured knockdown against classical micromechanics.

A composites journal will expect the computed law to be placed against the
standard estimates before it is offered as new. For spherical inclusions in an
isotropic matrix these are closed form:

  dilute (non-interacting Eshelby), Mori-Tanaka, and the Hashin-Shtrikman
  bounds. For a soft inclusion the relevant HS bound is the upper one.
"""
import numpy as np

# ---- phases (Section 3.4 of the paper)
E_m, nu_m = 9.4, 0.33
K_m = E_m / (3 * (1 - 2 * nu_m))
G_m = E_m / (2 * (1 + nu_m))

PHASES = {
    'brine (K=2.2 GPa, G=0.44 MPa)': (2.2, 0.44e-3),
    'gas void (K=0, G=0)':           (1e-9, 1e-9),
}


def E_of(K, G):
    return 9 * K * G / (3 * K + G)


def mori_tanaka(K_i, G_i, phi):
    """Mori-Tanaka estimate, spherical inclusions."""
    a = K_m + 4.0 / 3.0 * G_m
    f = G_m * (9 * K_m + 8 * G_m) / (6 * (K_m + 2 * G_m))
    K = K_m + phi * (K_i - K_m) * a / (a + (1 - phi) * (K_i - K_m))
    G = G_m + phi * (G_i - G_m) * (G_m + f) / ((G_m + f) + (1 - phi) * (G_i - G_m))
    return E_of(K, G)


def dilute(K_i, G_i, phi):
    """Non-interacting (dilute Eshelby) estimate, spherical inclusions."""
    a = K_m + 4.0 / 3.0 * G_m
    f = G_m * (9 * K_m + 8 * G_m) / (6 * (K_m + 2 * G_m))
    K = K_m + phi * (K_i - K_m) * a / (a + (K_i - K_m))
    G = G_m + phi * (G_i - G_m) * (G_m + f) / ((G_m + f) + (G_i - G_m))
    return E_of(K, G)


def hs_upper(K_i, G_i, phi):
    """Hashin-Shtrikman upper bound (stiff phase as comparison medium)."""
    c1, c2 = 1 - phi, phi            # 1 = matrix (stiff), 2 = inclusion (soft)
    K = K_m + c2 / (1.0 / (K_i - K_m) + c1 / (K_m + 4.0 / 3.0 * G_m))
    f = G_m * (9 * K_m + 8 * G_m) / (6 * (K_m + 2 * G_m))
    G = G_m + c2 / (1.0 / (G_i - G_m) + c1 / (G_m + f))
    return E_of(K, G)


print('matrix: E=%.2f GPa, nu=%.2f  ->  K=%.3f, G=%.3f GPa' % (E_m, nu_m, K_m, G_m))
print()
print('initial slope of E/E_m against phi (spherical inclusions)')
print('%-32s %9s %9s %9s   %s' % ('phase', 'dilute', 'Mori-T', 'HS upper', 'measured'))
print('-' * 82)

MEASURED = {'brine (K=2.2 GPa, G=0.44 MPa)': 1.68, 'gas void (K=0, G=0)': 1.64}

h = 1e-4
for name, (K_i, G_i) in PHASES.items():
    sl = {}
    for lab, fn in (('dilute', dilute), ('mt', mori_tanaka), ('hs', hs_upper)):
        sl[lab] = -(fn(K_i, G_i, h) / E_m - 1) / h
    print('%-32s %9.3f %9.3f %9.3f   %6.2f' %
          (name, sl['dilute'], sl['mt'], sl['hs'], MEASURED[name]))

print()
print('evaluated over the column range, E/E_m at phi = 0.02, 0.05, 0.10:')
print('%-14s %-10s %8s %8s %8s %8s' % ('phase', 'model', '0.02', '0.05', '0.10', '0.45'))
for name, (K_i, G_i) in PHASES.items():
    short = name.split()[0]
    for lab, fn in (('dilute', dilute), ('Mori-Tanaka', mori_tanaka), ('HS upper', hs_upper)):
        vals = [fn(K_i, G_i, p) / E_m for p in (0.02, 0.05, 0.10, 0.45)]
        print('%-14s %-10s %8.4f %8.4f %8.4f %8.4f' % (short, lab, *vals))
    m = MEASURED[name]
    vals = [1 - m * p for p in (0.02, 0.05, 0.10, 0.45)]
    print('%-14s %-10s %8.4f %8.4f %8.4f %8.4f' % (short, 'MEASURED', *vals))
    print()
