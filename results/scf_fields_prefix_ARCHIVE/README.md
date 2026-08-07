# Pre-fix stress fields (archived, do not use)

These .npz dumps were extracted from ODBs solved before the mesher fix, when
each ellipsoidal inclusion was built as its own bounding sphere. Their realised
soft-phase fractions are roughly twice the deck targets (e.g. WBL_BASE_s1 carries
phi_soft = 0.43 against a target of 0.21), which inflates every stress
concentration derived from them.

They are kept only so the discrepancy can be reproduced. The corrected
five-packing re-run lives in , and 
now holds copies of its s1 cells for the 3D field figure.
