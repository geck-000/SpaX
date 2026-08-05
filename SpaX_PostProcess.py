"""
SpaX_PostProcess.py
======================
Standalone post-processing for batch ODB extraction.
Uses odbAccess directly — NO CAE kernel needed.

Run with:
    abaqus python run_postprocess.py <csv_path> <odb_dir> [output_csv]
"""

from __future__ import print_function
import os
import sys
import csv
import numpy as np
import math
# NOTE: `from odbAccess import openOdb` is imported lazily inside the extract_*
# functions, not here. That keeps the merge mode (--merge) and the CLI dispatch
# importable under a plain `python3` (no Abaqus / no license needed) while the
# actual ODB reading still runs under `abaqus python`.


def mode_short(m):
    m = m.strip().lower()
    if 'uniaxial' in m and 'x' in m: return 'utx'
    if 'uniaxial' in m and 'y' in m: return 'uty'
    if 'uniaxial' in m and 'z' in m: return 'utz'
    if 'shear' in m and '12' in m: return 'ss12'
    if 'shear' in m and '13' in m: return 'ss13'
    if 'shear' in m and '23' in m: return 'ss23'
    if 'biaxial' in m: return 'bt'
    if 'confined' in m: return 'cc'
    return m.replace(' ', '')[:6]


def stress_component(m):
    m = m.strip().lower()
    if 'uniaxial' in m:
        if 'x' in m: return 'S11'
        if 'y' in m: return 'S22'
        if 'z' in m: return 'S33'
    if 'shear' in m:
        if '12' in m: return 'S12'
        if '13' in m: return 'S13'
        if '23' in m: return 'S23'
    return 'S11'


def stress_index(s_comp):
    idx_map = {'S11': 0, 'S22': 1, 'S33': 2, 'S12': 3, 'S13': 4, 'S23': 5}
    return idx_map.get(s_comp, 0)


def read_csv(csv_path):
    rows = []
    with open(csv_path, 'r') as f:
        for row in csv.DictReader(f):
            rows.append({k.strip(): (v.strip() if v else '') for k, v in row.items() if k})
    return rows


def _get_data(val):
    """Extract data array from an ODB field value."""
    try:
        return val.data
    except:
        return val.dataDouble


def _field_by_label(field, region):
    """Read an element field over `region` as (labels, data) NumPy arrays.

    Uses field.bulkDataBlocks -- ONE C-level call returning the whole field as
    arrays -- instead of iterating field.values element by element (each of
    which crosses the Abaqus C++/Python boundary). This is the dominant cost of
    post-processing, so the bulk read is ~10-100x faster. Falls back to .values
    if bulkDataBlocks is unavailable or empty (older releases / odd fields).

    Returns (labels: int ndarray (n,), data: float ndarray (n, ncomp)), with
    exactly ONE row per element: integration-point fields (C3D10/C3D10H have 4
    IPs/element) are collapsed to their per-element MEAN. This is essential for
    the whole-element EVOL volume-weighting the callers do — without it each of
    a quadratic element's 4 IP rows would be weighted by the FULL element volume
    (a 4x over-count of moments/volumes). For C3D4/C3D4H (1 IP) it is a no-op.
    """
    sub = field.getSubset(region=region)
    try:
        blocks = sub.bulkDataBlocks
    except Exception:
        blocks = None
    if blocks:
        labs, dats, ok = [], [], True
        for b in blocks:
            try:
                el, da = b.elementLabels, b.data
            except Exception:
                ok = False
                break
            if el is None or da is None:
                ok = False
                break
            d = np.asarray(da, dtype=float)
            if d.ndim == 1:
                d = d.reshape(-1, 1)
            labs.append(np.asarray(el).ravel())
            dats.append(d)
        if ok and labs:
            return _mean_per_element(np.concatenate(labs),
                                     np.concatenate(dats, axis=0))

    # Fallback: per-element .values (slow path).
    labs, dats = [], []
    for v in sub.values:
        labs.append(v.elementLabel)
        dats.append(np.asarray(_get_data(v), dtype=float).reshape(-1))
    if not labs:
        return np.array([], dtype=int), np.zeros((0, 1))
    return _mean_per_element(np.asarray(labs, dtype=int), np.vstack(dats))


def _mean_per_element(labels, data):
    """Collapse repeated-label rows (one per integration point) to one row per
    element, taking the mean over integration points. No-op when every label is
    already unique (1-IP elements)."""
    ulab, inv = np.unique(labels, return_inverse=True)
    if len(ulab) == len(labels):
        return labels, data
    sums = np.zeros((len(ulab), data.shape[1]), dtype=float)
    np.add.at(sums, inv, data)
    counts = np.bincount(inv).astype(float)
    return ulab, sums / counts[:, None]


def _aligned(values_map, labels):
    """Vectorized lookup: array of values_map[l] for l in labels, 0.0 if absent."""
    return np.fromiter((values_map.get(l, 0.0) for l in labels),
                       dtype=float, count=len(labels))


# Maps the driven stress component to the reference-point node set and the
# 0-based DOF on it that carries the imposed macro deformation. Normal modes
# drive RP-1/2/3 on their normal DOF; shear modes drive RP-4 on the shear DOF
# (ss12->2, ss13->1, ss23->2, matching SpaX_Standalone.shear_config).
RP_FOR_SCOMP = {
    'S11': ('RP-1', 0), 'S22': ('RP-2', 1), 'S33': ('RP-3', 2),
    'S12': ('RP-4', 1), 'S13': ('RP-4', 0), 'S23': ('RP-4', 1),
}


def _rp_region(odb, rp_name):
    """Assembly node set for a reference point, or None if absent."""
    sets = odb.rootAssembly.nodeSets
    for key in (rp_name, rp_name.upper(), rp_name.replace('-', '_').upper()):
        if key in sets.keys():
            return sets[key]
    return None


def _rp_disp(frame, rp_region, dof_idx):
    """U[dof_idx] of the single reference-point node in `frame`, or None.

    This is the actually-imposed macro deformation, independent of the load
    *Amplitude shape -- unlike the step-time proxy `frameValue * eng_strain`,
    which only equals it for a linear ramp.
    """
    if rp_region is None:
        return None
    try:
        uf = frame.fieldOutputs['U'].getSubset(region=rp_region)
        if len(uf.values) == 0:
            return None
        return float(_get_data(uf.values[0])[dof_idx])
    except Exception:
        return None


def extract_first_order(odb_path, s_comp, eng_strain, L):
    """
    Extract E or G and nu from an ODB using direct odbAccess.
    """
    from odbAccess import openOdb
    odb = openOdb(path=odb_path, readOnly=True)
    
    step = odb.steps[list(odb.steps.keys())[0]]
    s_idx = stress_index(s_comp)
    is_shear = 'S1' in s_comp and s_comp != 'S11' or s_comp in ('S12', 'S13', 'S23')
    
    V_RVE = L**3
    
    # Identify solid instances (exclude PBC_Surface)
    solid_instances = []
    for inst_name, inst in odb.rootAssembly.instances.items():
        if 'PBC' not in inst_name.upper() and inst_name.upper() != 'ASSEMBLY':
            solid_instances.append((inst_name, inst))
    
    # Reference point that carries the imposed macro deformation for this mode.
    rp_name, rp_dof = RP_FOR_SCOMP.get(s_comp, (None, None))
    rp_region = _rp_region(odb, rp_name) if rp_name else None

    stress_strain_data = []

    for frame_idx, frame in enumerate(step.frames):
        if frame_idx == 0:
            continue

        stress_field = frame.fieldOutputs['S']
        volume_field = frame.fieldOutputs['EVOL']
        
        total_stress_vol = 0.0
        total_vol = 0.0
        total_incl_vol = 0.0
        total_strain_axial_vol = 0.0
        total_strain_trans1_vol = 0.0
        total_strain_trans2_vol = 0.0
        
        # Check for strain field
        strain_key = None
        for key in ['LE', 'E']:
            if key in frame.fieldOutputs.keys():
                strain_key = key
                break
        
        for inst_name, inst in solid_instances:
            # Bulk-read S, EVOL, (LE/E) as arrays, then volume-weight with NumPy.
            # Same math as the old per-element loop, vectorized.
            s_lab, s_dat = _field_by_label(stress_field, inst)
            v_lab, v_dat = _field_by_label(volume_field, inst)
            if len(s_lab) == 0:
                continue
            labs = s_lab.tolist()
            vol_map = dict(zip(v_lab.tolist(), v_dat[:, 0].tolist()))
            w = _aligned(vol_map, labs)            # vol per stress-element, 0 if no vol

            total_stress_vol += float(np.dot(s_dat[:, s_idx], w))
            total_vol += float(w.sum())

            # Volume of the meshed soft phase, so the ACHIEVED inclusion
            # fraction can be reported. The void porosity below counts only
            # non-meshed volume (gas); brine is a meshed soft solid and is
            # therefore invisible to it, which makes void porosity useless for
            # auditing whether the packer reached its brine target.
            for _sname in ('SPHERE_ONLY', 'Sphere_Only'):
                try:
                    _set = inst.elementSets[_sname]
                except (KeyError, AttributeError):
                    continue
                try:
                    _l, _d = _field_by_label(volume_field, _set)
                    if len(_l):
                        total_incl_vol += float(_d[:, 0].sum())
                except Exception:
                    pass
                break

            if strain_key:
                e_lab, e_dat = _field_by_label(frame.fieldOutputs[strain_key], inst)
                # Align strain rows to the stress-element order; rows with no
                # strain (or no vol, via w=0) contribute nothing -- as before.
                erow = {int(l): i for i, l in enumerate(e_lab.tolist())}
                idx = np.fromiter((erow.get(l, -1) for l in labs),
                                  dtype=int, count=len(labs))
                estr = np.zeros((len(labs), e_dat.shape[1]))
                valid = idx >= 0
                if valid.any():
                    estr[valid] = e_dat[idx[valid]]
                total_strain_axial_vol += float(np.dot(estr[:, s_idx], w))
                if not is_shear:
                    t1, t2 = {'S11': (1, 2), 'S22': (0, 2), 'S33': (0, 1)}.get(
                        s_comp, (None, None))
                    if t1 is not None:
                        total_strain_trans1_vol += float(np.dot(estr[:, t1], w))
                        total_strain_trans2_vol += float(np.dot(estr[:, t2], w))
        
        # Hill-Mandel consistent:
        # Stress: sum(sigma*Ve) / V_RVE
        # Axial strain: RP-based (frameValue * eng_strain)
        # Transverse strain: element-averaged / V_solid (for Poisson only)
        sigma_macro = total_stress_vol / V_RVE
        # Macro strain = imposed RP displacement / L (Hill-Mandel macro strain
        # over the full RVE, voids included; amplitude-shape independent). Fall
        # back to the step-time proxy if RP output is unavailable in the ODB.
        rp_u = _rp_disp(frame, rp_region, rp_dof)
        eps_macro = (rp_u / L) if rp_u is not None else (frame.frameValue * eng_strain)
        
        if total_vol > 0 and not is_shear:
            eps_axial_solid = total_strain_axial_vol / total_vol
            eps_trans1 = total_strain_trans1_vol / total_vol
            eps_trans2 = total_strain_trans2_vol / total_vol
        else:
            eps_axial_solid = 0.0
            eps_trans1 = 0.0
            eps_trans2 = 0.0
        
        stress_strain_data.append({
            'sigma': sigma_macro,
            'eps': eps_macro,
            'eps_axial_solid': eps_axial_solid,
            'eps_trans1': eps_trans1,
            'eps_trans2': eps_trans2,
            # Meshed solid volume in this frame. Carried out of the loop so the
            # ACHIEVED porosity can be reported alongside the moduli: the packer
            # does not always reach the requested volume fraction, and a target
            # VoF echoed from the deck is therefore not evidence of what was
            # actually built. The second-order path already reports this; the
            # first-order path did not, which made high-VoF cases impossible to
            # audit from the results file alone.
            'v_solid': total_vol,
            'v_incl': total_incl_vol,
        })

    odb.close()

    results = {}
    if stress_strain_data and V_RVE > 0:
        # Small-strain volume change is negligible, but average over frames
        # rather than trusting any single one.
        v_solid = sum(d['v_solid'] for d in stress_strain_data) / float(len(stress_strain_data))
        v_incl = sum(d['v_incl'] for d in stress_strain_data) / float(len(stress_strain_data))
        results['V_solid'] = v_solid
        # Non-meshed volume only. For these decks that is the gas, since brine
        # is a meshed soft phase and therefore counts as solid here.
        results['porosity'] = 1.0 - v_solid / V_RVE
        # Achieved soft-phase (brine) fraction, and the two together -- the
        # quantity to compare against the requested VoF when auditing whether
        # the packer met its target.
        results['phi_inclusion'] = v_incl / V_RVE
        results['phi_soft_total'] = (v_incl / V_RVE) + (1.0 - v_solid / V_RVE)
    if len(stress_strain_data) >= 2:
        eps_arr = np.array([d['eps'] for d in stress_strain_data])
        sig_arr = np.array([d['sigma'] for d in stress_strain_data])
        
        # Linear regression over linear region (matching kernel's polyfit)
        # Use points from 10% to 40% of max strain (same as kernel)
        max_eps = max(abs(eps_arr))
        if max_eps > 0:
            linear_mask = (np.abs(eps_arr) >= 0.1 * max_eps) & (np.abs(eps_arr) <= 0.4 * max_eps)
            if np.sum(linear_mask) >= 2:
                coeffs = np.polyfit(eps_arr[linear_mask], sig_arr[linear_mask], 1)
                modulus = coeffs[0]
            else:
                # Fallback: use all points
                coeffs = np.polyfit(eps_arr, sig_arr, 1)
                modulus = coeffs[0]
        else:
            modulus = 0
        
        if is_shear:
            results['G_eff'] = modulus
        else:
            results['E_eff'] = modulus
            
            # Poisson's ratio from linear regression (matching kernel)
            if np.sum(linear_mask) >= 2:
                eps_lin = eps_arr[linear_mask]
                axial_solid_lin = np.array([stress_strain_data[j]['eps_axial_solid'] 
                                       for j in range(len(stress_strain_data)) if linear_mask[j]])
                trans1_lin = np.array([stress_strain_data[j]['eps_trans1'] 
                                       for j in range(len(stress_strain_data)) if linear_mask[j]])
                trans2_lin = np.array([stress_strain_data[j]['eps_trans2'] 
                                       for j in range(len(stress_strain_data)) if linear_mask[j]])
                nu1 = -np.polyfit(axial_solid_lin, trans1_lin, 1)[0]
                nu2 = -np.polyfit(axial_solid_lin, trans2_lin, 1)[0]
                results['nu_eff'] = (nu1 + nu2) / 2.0
            else:
                # Fallback: last frame
                if abs(eps_arr[-1]) > 1e-30:
                    nu1 = -stress_strain_data[-1]['eps_trans1'] / stress_strain_data[-1]['eps_axial_solid']
                    nu2 = -stress_strain_data[-1]['eps_trans2'] / stress_strain_data[-1]['eps_axial_solid']
                    results['nu_eff'] = (nu1 + nu2) / 2.0
    
    return results


def extract_second_order(odb_path, L, Kappa, Bending_Plane):
    """
    Extract D_RVE and E_bending from a bending ODB using direct odbAccess.
    """
    from odbAccess import openOdb
    odb = openOdb(path=odb_path, readOnly=True)
    
    # Find bending step
    step_name = None
    for sname in odb.steps.keys():
        if 'Bending' in sname or 'BENDING' in sname:
            step_name = sname
            break
    if step_name is None:
        step_name = list(odb.steps.keys())[0]
    step = odb.steps[step_name]
    
    V_RVE = L**3
    
    # Stress/position indices by bending plane
    if Bending_Plane == 'xz':
        s_idx = 0; coord_idx = 2; coord_bar = L / 2.0
    elif Bending_Plane == 'yz':
        s_idx = 1; coord_idx = 2; coord_bar = L / 2.0
    elif Bending_Plane == 'xy':
        s_idx = 0; coord_idx = 1; coord_bar = L / 2.0
    
    # Solid instances
    solid_instances = []
    for inst_name, inst in odb.rootAssembly.instances.items():
        if 'PBC' not in inst_name.upper() and inst_name.upper() != 'ASSEMBLY':
            solid_instances.append((inst_name, inst))
    
    # Build centroid lookup
    elem_centroid = {}
    for inst_name, inst in solid_instances:
        node_coords = {}
        for node in inst.nodes:
            node_coords[node.label] = node.coordinates
        for elem in inst.elements:
            conn = elem.connectivity
            n = len(conn)
            cx, cy, cz = 0.0, 0.0, 0.0
            for nid in conn:
                c = node_coords[nid]
                cx += c[0]; cy += c[1]; cz += c[2]
            elem_centroid[(inst_name, elem.label)] = (cx/n, cy/n, cz/n)
    
    # Process last frame. A bending analysis that failed to converge writes an
    # ODB with an empty step (0 frames); guard against it so the caller records
    # this cleanly instead of dying on an out-of-range frame index.
    if len(step.frames) == 0:
        odb.close()
        raise RuntimeError(
            "bending step '{}' has 0 frames (analysis did not converge)".format(step_name))
    last_frame = step.frames[-1]
    stress_field = last_frame.fieldOutputs['S']
    volume_field = last_frame.fieldOutputs['EVOL']
    
    sigma_vol = np.zeros(6)
    moment_vol = 0.0
    total_vol = 0.0
    
    for inst_name, inst in solid_instances:
        # Bulk-read S and EVOL as arrays, then volume-weight with NumPy.
        s_lab, s_dat = _field_by_label(stress_field, inst)
        v_lab, v_dat = _field_by_label(volume_field, inst)
        if len(s_lab) == 0:
            continue
        labs = s_lab.tolist()
        vol_map = dict(zip(v_lab.tolist(), v_dat[:, 0].tolist()))

        # Weight is the element volume, but ZERO for elements missing a volume
        # OR a centroid (the old loop `continue`d on either).
        cent_z = {el: c[coord_idx] for (iname, el), c in elem_centroid.items()
                  if iname == inst_name}
        w = _aligned(vol_map, labs)
        in_cent = np.fromiter((1.0 if l in cent_z else 0.0 for l in labs),
                              dtype=float, count=len(labs))
        w = w * in_cent
        z_rel = _aligned(cent_z, labs) - coord_bar

        ncomp = min(6, s_dat.shape[1])
        sigma_vol[:ncomp] += (s_dat[:, :ncomp] * w[:, None]).sum(axis=0)
        moment_vol += float(np.dot(s_dat[:, s_idx] * z_rel, w))
        total_vol += float(w.sum())
    
    # Imposed curvature: read RP_K's prescribed DOF directly (amplitude-shape
    # independent), falling back to the step-time proxy if RP_K is unavailable.
    rpk_u = _rp_disp(last_frame, _rp_region(odb, 'RP_K'), 0)
    kappa_actual = rpk_u if (rpk_u is not None and abs(rpk_u) > 1e-30) \
        else Kappa * last_frame.frameValue
    M_over_V = moment_vol / L
    
    D_rve = abs(M_over_V / kappa_actual) if abs(kappa_actual) > 1e-30 else 0
    E_bending = 12.0 * D_rve / L**4
    sigma_bar = sigma_vol / V_RVE
    porosity = 1.0 - total_vol / V_RVE
    N_membrane = sigma_vol[s_idx] / L
    B_coupling = N_membrane / kappa_actual if abs(kappa_actual) > 1e-30 else 0
    
    odb.close()
    
    return {
        'D_rve': D_rve,
        'E_bending': E_bending,
        'porosity': porosity,
        'sigma_bar': sigma_bar,
        'N_membrane': N_membrane,
        'B_coupling': B_coupling,
    }


def compute_length_scale(D_rve, E_eff, G_eff, nu_eff, L):
    """
    Compute MCST length scale from D_RVE and first-order properties.
    Matches kernel's extract_length_scale exactly.
    """
    plate_factor = 1.0 / (1.0 - nu_eff**2)
    E_plate = E_eff * plate_factor
    I_cs = L**4 / 12.0
    D_classical = E_plate * I_cs
    
    l_squared = (D_rve - D_classical) / (G_eff * L**2) if G_eff > 0 else 0
    l = (abs(l_squared))**0.5 if l_squared > 0 else 0.0
    
    E_bending_plate = D_rve / I_cs if I_cs > 0 else 0.0
    E_bending_material = E_bending_plate * (1.0 - nu_eff**2)
    
    return {
        'D_classical': D_classical,
        'l_squared': l_squared,
        'l': l,
        'D_ratio': D_rve / D_classical if D_classical > 0 else 0,
        'E_bending_material': E_bending_material,
        'E_bending_plate': E_bending_plate,
    }


def run(csv_path, odb_dir, output_csv='postprocess_results.csv', only_index=None):
    """Main post-processing pipeline.

    only_index : 1-based row index into the CSV. When set, process ONLY that one
    RVE and write a 1-row CSV. This is what lets stage 3 run as a SLURM job array
    (one task per RVE) for parallel post-processing; a separate merge step
    (merge_partials) then unions the per-task CSVs into the final results.csv.
    When None, every RVE is processed in one pass (original behaviour).
    """
    params_list = read_csv(csv_path)

    if only_index is not None:
        if only_index < 1 or only_index > len(params_list):
            raise IndexError("only_index {} out of range 1..{}".format(
                only_index, len(params_list)))
        params_list = [params_list[only_index - 1]]

    print("\n" + "=" * 70)
    print("SPAX POST-PROCESSING PIPELINE")
    print("=" * 70)
    print("  CSV: {}".format(csv_path))
    print("  ODB dir: {}".format(odb_dir))
    if only_index is not None:
        print("  RVE: #{} ({})".format(only_index, params_list[0].get('run_id', '?')))
    else:
        print("  RVEs: {}".format(len(params_list)))
    print("  Output: {}".format(output_csv))

    all_results = []
    
    for i, params in enumerate(params_list):
        run_id = params.get('run_id', 'rve')
        L = float(params['L'])
        Mode = params.get('Mode', 'Uniaxial Tension X')
        Mode2 = params.get('Mode2', '')
        Disp = float(params.get('Disp', L * 0.01))
        Disp2 = float(params.get('Disp2', Disp)) if params.get('Disp2', '') else Disp
        Kappa = float(params.get('Kappa', 0.0))
        Bending_Plane = params.get('Bending_Plane', 'xz')
        
        print("\n" + "-" * 70)
        print("[{}/{}] {}".format(i + 1, len(params_list), run_id))
        print("-" * 70)
        
        row = {'run_id': run_id}
        for key in ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix',
                     'VoF_sphere', 'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere',
                     'E_sphere_inclusion', 'sphericity_avg', 'PBC_Method',
                     'Bending_PBC_Type', 'Growth_Direction', 'generate_channels',
                     'channel_vof_target']:
            row[key] = params.get(key, '')
        
        # ---- First-order: directional moduli (ANISOTROPY) ----
        # Probe every canonical first-order ODB present and extract whichever
        # were solved. An isotropic run solves only Mode(+Mode2); a full-tensor
        # run solves all six and yields E_x/E_y/E_z and G_xy/G_xz/G_yz, so the
        # RVE's ANISOTROPY is visible -- e.g. the transverse isotropy of the
        # vertical (Z) brine-channel network, where E_z != E_x ~ E_y. The
        # backward-compatible scalars E_eff/nu_eff/G_eff (consumed by the MCST
        # length-scale block below and by existing CSV readers) are then taken
        # from the run's declared Mode/Mode2, falling back to the first present.
        eng = (Disp / L) if L > 0 else 0.01
        UNI = [('utx', 'S11', 'E_x', 'nu_x'), ('uty', 'S22', 'E_y', 'nu_y'),
               ('utz', 'S33', 'E_z', 'nu_z')]
        SHR = [('ss12', 'S12', 'G_xy'), ('ss13', 'S13', 'G_xz'),
               ('ss23', 'S23', 'G_yz')]

        def _fo(short, scomp):
            p = os.path.join(odb_dir, 'Job-{}-{}.odb'.format(run_id, short))
            if not os.path.isfile(p):
                return None
            try:
                return extract_first_order(p, scomp, eng, L)
            except Exception as e:
                print("    ERROR {}: {}".format(short, e))
                return 'ERROR'

        for short, scomp, Ek, nuk in UNI:
            r = _fo(short, scomp)
            if r == 'ERROR':
                row[Ek] = 'ERROR'
            elif r:
                row[Ek] = r.get('E_eff', ''); row[nuk] = r.get('nu_eff', '')
                # Achieved porosity, from the meshed solid volume rather than
                # the requested VoF. Recorded from the first uniaxial ODB that
                # yields it (all load cases share the geometry), and only when
                # the bending path has not already supplied it.
                if 'porosity' in r and row.get('porosity', '') in ('', None):
                    row['porosity'] = r['porosity']
                    row['V_solid'] = r.get('V_solid', '')
                    row['phi_inclusion'] = r.get('phi_inclusion', '')
                    row['phi_soft_total'] = r.get('phi_soft_total', '')
                print("    {} -> E={:.4e}  nu={}".format(short, r.get('E_eff', 0),
                                                         r.get('nu_eff', '')))
        for short, scomp, Gk in SHR:
            r = _fo(short, scomp)
            if r == 'ERROR':
                row[Gk] = 'ERROR'
            elif r:
                row[Gk] = r.get('G_eff', r.get('E_eff', ''))
                try:
                    print("    {} -> G={:.4e}".format(short, float(row[Gk])))
                except (TypeError, ValueError):
                    pass

        def _pick(keymap, default_first):
            """Value for the declared mode's column, else the first present."""
            k = keymap.get(default_first)
            if k and row.get(k, '') not in ('', None):
                return row[k]
            for kk in keymap.values():
                if row.get(kk, '') not in ('', None):
                    return row[kk]
            return 'MISSING'

        row['E_eff'] = _pick({'utx': 'E_x', 'uty': 'E_y', 'utz': 'E_z'},
                             mode_short(Mode))
        row['nu_eff'] = _pick({'utx': 'nu_x', 'uty': 'nu_y', 'utz': 'nu_z'},
                              mode_short(Mode))
        row['G_eff'] = _pick({'ss12': 'G_xy', 'ss13': 'G_xz', 'ss23': 'G_yz'},
                             mode_short(Mode2) if Mode2 else 'ss13')

        # ---- Anisotropy summary ----
        def _num(k):
            v = row.get(k, '')
            try:
                return float(v) if v not in ('', 'ERROR', 'MISSING') else None
            except (TypeError, ValueError):
                return None
        Es = [v for v in (_num('E_x'), _num('E_y'), _num('E_z')) if v and v > 0]
        if len(Es) >= 2:
            row['E_anisotropy'] = max(Es) / min(Es)        # 1.0 == isotropic
            Ex, Ey, Ez = _num('E_x'), _num('E_y'), _num('E_z')
            if Ex and Ey and Ez:
                # transverse-isotropy index about Z (vertical channels): axial Ez
                # vs in-plane mean. >1 stiffer along the channels, <1 softer.
                row['E_z_over_xy'] = Ez / (0.5 * (Ex + Ey))
            print("    Anisotropy: E_max/E_min = {:.3f}{}".format(
                row['E_anisotropy'],
                "  E_z/E_xy = {:.3f}".format(row['E_z_over_xy'])
                if 'E_z_over_xy' in row else ""))
        
        # ---- Second-order: Bending ----
        odb_ben = os.path.join(odb_dir, 'Job-{}-ben.odb'.format(run_id))
        if Kappa > 0 and os.path.isfile(odb_ben):
            print("  Extracting bending (D_RVE, E_bending)...")
            try:
                r3 = extract_second_order(odb_ben, L, Kappa, Bending_Plane)
                row['D_rve'] = r3.get('D_rve', '')
                row['E_bending'] = r3.get('E_bending', '')
                row['porosity'] = r3.get('porosity', '')
                row['N_membrane'] = r3.get('N_membrane', '')
                row['B_coupling'] = r3.get('B_coupling', '')
                print("    D_RVE = {:.4e}, E_bend = {:.4e}, porosity = {:.2%}".format(
                    r3.get('D_rve', 0), r3.get('E_bending', 0), r3.get('porosity', 0)))
            except Exception as e:
                print("    ERROR (2nd order): {}".format(e))
                row['D_rve'] = 'ERROR'
        elif Kappa > 0:
            print("  [SKIP] {} not found".format(os.path.basename(odb_ben)))
            row['D_rve'] = 'MISSING'
        
        # ---- MCST Length Scale ----
        E_eff = float(row.get('E_eff', 0) or 0) if row.get('E_eff', '') not in ('ERROR', 'MISSING', '') else 0
        G_eff = float(row.get('G_eff', 0) or 0) if row.get('G_eff', '') not in ('ERROR', 'MISSING', '') else 0
        D_rve = float(row.get('D_rve', 0) or 0) if row.get('D_rve', '') not in ('ERROR', 'MISSING', '') else 0
        nu_eff = float(row.get('nu_eff', 0.3) or 0.3) if row.get('nu_eff', '') not in ('ERROR', 'MISSING', '') else 0.3
        
        if D_rve > 0 and E_eff > 0 and G_eff > 0:
            print("  Computing MCST length scale...")
            try:
                rl = compute_length_scale(D_rve, E_eff, G_eff, nu_eff, L)
                row['D_classical'] = rl['D_classical']
                row['l_squared'] = rl['l_squared']
                row['l'] = rl['l']
                row['D_ratio'] = rl['D_ratio']
                row['E_bending_material'] = rl['E_bending_material']
                print("    D_class = {:.4e}, D/D_class = {:.4f}, l = {:.4e}, l/L = {:.4f}".format(
                    rl['D_classical'], rl['D_ratio'], rl['l'], rl['l'] / L if L > 0 else 0))
            except Exception as e:
                print("    ERROR (MCST): {}".format(e))
                row['l'] = 'ERROR'
        elif Kappa > 0 and (E_eff <= 0 or G_eff <= 0):
            print("  [SKIP] MCST: missing E_eff or G_eff")
        
        all_results.append(row)
    
    # Write consolidated CSV
    if all_results:
        all_keys = []
        for r in all_results:
            for k in r.keys():
                if k not in all_keys:
                    all_keys.append(k)
        
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_keys)
            writer.writeheader()
            writer.writerows(all_results)
        
        print("\n" + "=" * 70)
        print("POST-PROCESSING COMPLETE")
        print("  {} RVEs processed".format(len(all_results)))
        print("  Results: {}".format(output_csv))
        print("=" * 70)
    
    return all_results


def print_summary_table(results):
    print("\n" + "=" * 120)
    print("{:<15} {:>12} {:>12} {:>8} {:>12} {:>12} {:>12} {:>8} {:>8}".format(
        'run_id', 'E_eff', 'G_eff', 'nu_eff', 'D_RVE', 'D/D_class', 'l', 'l/L', 'porosity'))
    print("-" * 120)
    
    for r in results:
        run_id = r.get('run_id', '?')
        L = float(r.get('L', 0.2) or 0.2)
        
        def fmt(val, w=12):
            if val in ('', 'ERROR', 'MISSING', 'MISSING_PROPS', None):
                return '{:>{w}}'.format('-', w=w)
            try:
                v = float(val)
                if abs(v) > 1e6: return '{:>{w}.4e}'.format(v, w=w)
                elif abs(v) > 0.001: return '{:>{w}.4f}'.format(v, w=w)
                else: return '{:>{w}.4e}'.format(v, w=w)
            except:
                return '{:>{w}}'.format(str(val)[:w], w=w)
        
        l_val = r.get('l', '')
        try:
            l_over_L = '{:.4f}'.format(float(l_val) / L)
        except:
            l_over_L = '-'
        
        print("{:<15} {} {} {:>8} {} {} {} {:>8} {}".format(
            run_id[:15],
            fmt(r.get('E_eff', '')),
            fmt(r.get('G_eff', '')),
            fmt(r.get('nu_eff', ''), 8),
            fmt(r.get('D_rve', '')),
            fmt(r.get('D_ratio', '')),
            fmt(r.get('l', '')),
            l_over_L,
            fmt(r.get('porosity', ''), 8)))
    
    print("=" * 120)


def merge_partials(parts_dir, output_csv):
    """Union the per-RVE partial CSVs written by the post-processing job array
    into one results.csv. Columns vary by RVE (e.g. only bending rows have
    D_rve), so take the union of all columns (first-seen order) and fill blanks.
    Rows are ordered by the row index encoded in each partial's filename so the
    output matches the input CSV order. Pure-Python: run with plain python3."""
    import glob
    import os
    import re

    files = glob.glob(os.path.join(parts_dir, 'row_*.csv'))

    def _idx(f):
        m = re.search(r'row_(\d+)\.csv$', os.path.basename(f))
        return int(m.group(1)) if m else (1 << 30)

    files.sort(key=_idx)

    rows = []
    keys = []
    for f in files:
        with open(f) as fh:
            for r in csv.DictReader(fh):
                rows.append(r)
                for k in r.keys():
                    if k not in keys:
                        keys.append(k)

    if not rows:
        print("merge: no partial rows found in {}".format(parts_dir))
        return []

    with open(output_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})

    print("merged {} rows from {} partial file(s) -> {}".format(
        len(rows), len(files), output_csv))
    return rows




# =====================================================================
# PRINCIPAL / FULL-TENSOR / ELASTICITY-TENSOR EXTRACTION  (abaqus python).
# Folded in from the former SpaX_PrincipalExtract.py. The 6x6 elasticity
# tensor C_ij is the COMPLETE anisotropy characterisation of an RVE.
#   abaqus python SpaX_PostProcess.py principal   <odb> <out> [L] [run_id]
#   abaqus python SpaX_PostProcess.py elasticity  <odb_dir> <out> <L> <run_id>
# =====================================================================
def _openOdb(*a, **k):
    from odbAccess import openOdb
    return openOdb(*a, **k)


def _eigen_symmetric_3x3(s11, s22, s33, s12, s13, s23):
    """
    Compute eigenvalues (principal values) of a 3x3 symmetric tensor.
    Uses numpy for robustness.
    Returns (p1, p2, p3) sorted: p1 >= p2 >= p3.
    """
    mat = np.array([
        [s11, s12, s13],
        [s12, s22, s23],
        [s13, s23, s33]
    ])
    eigvals = np.linalg.eigvalsh(mat)
    return float(eigvals[2]), float(eigvals[1]), float(eigvals[0])



def extract_principals(Odb_Path, Output_Path, L=0.0, run_id='RVE',
                        last_frame_only=False, subtract_mean=False):
    """
    Extract principal stresses, strains, strain energy density,
    and full volume-averaged stress/strain tensors from an RVE ODB.
    
    Parameters
    ----------
    Odb_Path : str
        Full path to the ODB file.
    Output_Path : str
        Directory to write CSV results.
    L : float
        RVE side length. If > 0, volume-average uses V_RVE = L^3
        (correct for porous). If 0, uses sum of element volumes.
    run_id : str
        Identifier for output filenames.
    last_frame_only : bool
        If True, only process the last frame (faster).
    subtract_mean : bool
        If True, subtract the mean principal stress and strain from
        the volume-averaged values. Use for bending load cases to
        remove the coordinate-shift artefact from non-centred PBCs
        (origin at corner instead of midplane).
    
    Returns
    -------
    dict with volume-averaged principal quantities and full tensors
    at last frame.
    """
    
    if not os.path.isfile(Odb_Path):
        raise FileNotFoundError("ODB not found: " + Odb_Path)
    
    odb = _openOdb(path=Odb_Path, readOnly=True)
    
    # Pick the step with the most frames (handles suppressed steps)
    step_name = None
    max_frames = -1
    for sn in odb.steps.keys():
        nf = len(odb.steps[sn].frames)
        if nf > max_frames:
            max_frames = nf
            step_name = sn
    
    if step_name is None or max_frames <= 0:
        print("  ERROR: ODB has no steps with frames: {}".format(Odb_Path))
        print("  Available steps: {}".format(list(odb.steps.keys())))
        odb.close()
        return None
    
    step = odb.steps[step_name]
    
    print("\n" + "=" * 70)
    print("PRINCIPAL STRESS/STRAIN/ENERGY EXTRACTION")
    print("=" * 70)
    print("  ODB: {}".format(Odb_Path))
    print("  Step: {} ({} frames)".format(step_name, len(step.frames)))
    
    # Identify solid instances (exclude PBC_Surface and ASSEMBLY pseudo-instance)
    solid_instances = []
    for inst_name in odb.rootAssembly.instances.keys():
        inst_upper = inst_name.upper()
        if 'PBC' in inst_upper:
            continue
        if inst_upper == 'ASSEMBLY':
            continue
        solid_instances.append(inst_name)
    print("  Solid instances: {}".format(', '.join(solid_instances)))
    
    # Check available fields
    last_frame = step.frames[-1]
    available = list(last_frame.fieldOutputs.keys())
    has_stress = 'S' in available
    has_evol = 'EVOL' in available
    
    strain_key = 'LE' if 'LE' in available else ('E' if 'E' in available else None)
    
    print("  Available: S={}, strain={}, EVOL={}".format(
        has_stress, strain_key, has_evol))
    
    V_RVE = L**3 if L > 0 else 0.0
    
    # Determine which frames to process
    if last_frame_only:
        frame_indices = [len(step.frames) - 1]
    else:
        frame_indices = range(len(step.frames))
    
    frame_results = []
    
    for fi in frame_indices:
        frame = step.frames[fi]
        frame_time = frame.frameValue
        
        stress_field = frame.fieldOutputs['S'] if has_stress else None
        strain_field = frame.fieldOutputs[strain_key] if strain_key else None
        vol_field = frame.fieldOutputs['EVOL'] if has_evol else None
        
        # Principal accumulators
        sum_sp1_v = 0.0; sum_sp2_v = 0.0; sum_sp3_v = 0.0
        sum_ep1_v = 0.0; sum_ep2_v = 0.0; sum_ep3_v = 0.0
        sum_vm_v = 0.0; sum_press_v = 0.0; sum_tresca_v = 0.0
        sum_se_v = 0.0
        
        # Full tensor accumulators (Voigt: 11, 22, 33, 12, 13, 23)
        sum_s_v = np.zeros(6)
        sum_e_v = np.zeros(6)
        
        sum_vol = 0.0
        
        all_sp1 = []; all_sp2 = []; all_sp3 = []
        all_vm = []; all_press = []
        
        for inst_name in solid_instances:
            s_vals = []
            e_vals = []
            v_vals = []
            
            if stress_field:
                try:
                    s_vals = stress_field.getSubset(
                        region=odb.rootAssembly.instances[inst_name]).values
                except:
                    pass
            
            if strain_field:
                try:
                    e_vals = strain_field.getSubset(
                        region=odb.rootAssembly.instances[inst_name]).values
                except:
                    pass
            
            if vol_field:
                try:
                    v_vals = vol_field.getSubset(
                        region=odb.rootAssembly.instances[inst_name]).values
                except:
                    pass
            
            n_elem = len(s_vals)
            if fi == frame_indices[-1]:
                inst_vol = sum([abs(_get_data(v_vals[i])) for i in range(min(n_elem, len(v_vals)))]) if v_vals else 0.0
                print("    Instance {}: {} elements, V={:.6e}".format(inst_name, n_elem, inst_vol))
            
            for i in range(n_elem):
                ve = abs(_get_data(v_vals[i])) if i < len(v_vals) else 1.0
                
                if i < len(s_vals):
                    sd = _get_data(s_vals[i])
                    s11 = float(sd[0]); s22 = float(sd[1]); s33 = float(sd[2])
                    s12 = float(sd[3]); s13 = float(sd[4]); s23 = float(sd[5])
                    
                    # Full stress tensor
                    sum_s_v[0] += s11 * ve
                    sum_s_v[1] += s22 * ve
                    sum_s_v[2] += s33 * ve
                    sum_s_v[3] += s12 * ve
                    sum_s_v[4] += s13 * ve
                    sum_s_v[5] += s23 * ve
                    
                    sp1, sp2, sp3 = _eigen_symmetric_3x3(
                        s11, s22, s33, s12, s13, s23)
                    vm = math.sqrt(0.5 * ((sp1-sp2)**2 + (sp2-sp3)**2 + (sp3-sp1)**2))
                    press = -(s11 + s22 + s33) / 3.0
                    tresca = sp1 - sp3
                    
                    sum_sp1_v += sp1 * ve
                    sum_sp2_v += sp2 * ve
                    sum_sp3_v += sp3 * ve
                    sum_vm_v += vm * ve
                    sum_press_v += press * ve
                    sum_tresca_v += tresca * ve
                    
                    all_sp1.append(sp1)
                    all_sp2.append(sp2)
                    all_sp3.append(sp3)
                    all_vm.append(vm)
                    all_press.append(press)
                
                if i < len(e_vals):
                    ed = _get_data(e_vals[i])
                    e11 = float(ed[0]); e22 = float(ed[1]); e33 = float(ed[2])
                    e12 = float(ed[3]); e13 = float(ed[4]); e23 = float(ed[5])
                    
                    # Full strain tensor
                    sum_e_v[0] += e11 * ve
                    sum_e_v[1] += e22 * ve
                    sum_e_v[2] += e33 * ve
                    sum_e_v[3] += e12 * ve
                    sum_e_v[4] += e13 * ve
                    sum_e_v[5] += e23 * ve
                    
                    ep1, ep2, ep3 = _eigen_symmetric_3x3(
                        e11, e22, e33, e12, e13, e23)
                    sum_ep1_v += ep1 * ve
                    sum_ep2_v += ep2 * ve
                    sum_ep3_v += ep3 * ve
                
                if i < len(s_vals) and i < len(e_vals):
                    sd2 = _get_data(s_vals[i])
                    ed2 = _get_data(e_vals[i])
                    sed = 0.5 * (float(sd2[0])*float(ed2[0]) + 
                                 float(sd2[1])*float(ed2[1]) + 
                                 float(sd2[2])*float(ed2[2]) +
                                 2.0*float(sd2[3])*float(ed2[3]) + 
                                 2.0*float(sd2[4])*float(ed2[4]) + 
                                 2.0*float(sd2[5])*float(ed2[5]))
                    sum_se_v += sed * ve
                
                sum_vol += ve
        
        # Volume averaging
        V_avg = V_RVE if V_RVE > 0 else sum_vol
        
        if fi == frame_indices[0] or fi == frame_indices[-1]:
            porosity = 1.0 - sum_vol / (L**3) if L > 0 else 0.0
            print("    Frame {}: V_solid={:.6e}, V_RVE={:.6e}, porosity={:.1%}, using V={}".format(
                fi, sum_vol, V_avg, porosity,
                'L^3' if (V_RVE > 0) else 'V_solid'))
        
        if V_avg > 0:
            avg_sp1 = sum_sp1_v / V_avg
            avg_sp2 = sum_sp2_v / V_avg
            avg_sp3 = sum_sp3_v / V_avg
            avg_ep1 = sum_ep1_v / V_avg
            avg_ep2 = sum_ep2_v / V_avg
            avg_ep3 = sum_ep3_v / V_avg
            avg_vm = sum_vm_v / V_avg
            avg_press = sum_press_v / V_avg
            avg_tresca = sum_tresca_v / V_avg
            avg_se = sum_se_v / V_avg
            total_se = sum_se_v
            avg_stress = sum_s_v / V_avg
            avg_strain = sum_e_v / V_avg
        else:
            avg_sp1 = avg_sp2 = avg_sp3 = 0.0
            avg_ep1 = avg_ep2 = avg_ep3 = 0.0
            avg_vm = avg_press = avg_tresca = 0.0
            avg_se = total_se = 0.0
            avg_stress = np.zeros(6)
            avg_strain = np.zeros(6)
        
        # Bending correction
        if subtract_mean:
            mean_s = (avg_sp1 + avg_sp2 + avg_sp3) / 3.0
            mean_e = (avg_ep1 + avg_ep2 + avg_ep3) / 3.0
            
            if fi == frame_indices[-1]:
                print("    Bending correction: subtracting mean stress = {:.4e}".format(mean_s))
                print("    Bending correction: subtracting mean strain = {:.6e}".format(mean_e))
            
            avg_sp1 -= mean_s
            avg_sp2 -= mean_s
            avg_sp3 -= mean_s
            avg_ep1 -= mean_e
            avg_ep2 -= mean_e
            avg_ep3 -= mean_e
            
            avg_stress[0] -= mean_s
            avg_stress[1] -= mean_s
            avg_stress[2] -= mean_s
            avg_strain[0] -= mean_e
            avg_strain[1] -= mean_e
            avg_strain[2] -= mean_e
            
            avg_press = -(avg_sp1 + avg_sp2 + avg_sp3) / 3.0
            avg_vm = math.sqrt(0.5 * ((avg_sp1-avg_sp2)**2 + 
                                       (avg_sp2-avg_sp3)**2 + 
                                       (avg_sp3-avg_sp1)**2))
            avg_tresca = avg_sp1 - avg_sp3
            avg_se -= abs(mean_s * mean_e)
            total_se = avg_se * V_avg
        
        result = {
            'frame': fi,
            'time': frame_time,
            'V_solid': sum_vol,
            'V_RVE': V_avg,
            # Full volume-averaged stress tensor (Voigt)
            'S11': float(avg_stress[0]),
            'S22': float(avg_stress[1]),
            'S33': float(avg_stress[2]),
            'S12': float(avg_stress[3]),
            'S13': float(avg_stress[4]),
            'S23': float(avg_stress[5]),
            # Full volume-averaged strain tensor (Voigt)
            'E11': float(avg_strain[0]),
            'E22': float(avg_strain[1]),
            'E33': float(avg_strain[2]),
            'E12': float(avg_strain[3]),
            'E13': float(avg_strain[4]),
            'E23': float(avg_strain[5]),
            # Volume-averaged principal stresses
            'sigma_1': avg_sp1,
            'sigma_2': avg_sp2,
            'sigma_3': avg_sp3,
            'sigma_vm': avg_vm,
            'pressure': avg_press,
            'tresca': avg_tresca,
            # Volume-averaged principal strains
            'epsilon_1': avg_ep1,
            'epsilon_2': avg_ep2,
            'epsilon_3': avg_ep3,
            # Strain energy
            'W_density': avg_se,
            'W_total': total_se,
            # Statistics
            'sp1_max': max(all_sp1) if all_sp1 else 0.0,
            'sp1_min': min(all_sp1) if all_sp1 else 0.0,
            'sp3_max': max(all_sp3) if all_sp3 else 0.0,
            'sp3_min': min(all_sp3) if all_sp3 else 0.0,
            'vm_max': max(all_vm) if all_vm else 0.0,
            'vm_mean': np.mean(all_vm) if all_vm else 0.0,
            'vm_std': np.std(all_vm) if all_vm else 0.0,
            'press_max': max(all_press) if all_press else 0.0,
            'press_min': min(all_press) if all_press else 0.0,
        }
        frame_results.append(result)
    
    odb.close()
    
    # ================================================================
    # Print summary
    # ================================================================
    
    last = frame_results[-1]
    
    print("\n" + "=" * 70)
    print("RESULTS (last frame, t={:.4f})".format(last['time']))
    print("=" * 70)
    
    print("\n  Volume-averaged STRESS TENSOR (Voigt):")
    print("    S11={:.6e}  S22={:.6e}  S33={:.6e}".format(
        last['S11'], last['S22'], last['S33']))
    print("    S12={:.6e}  S13={:.6e}  S23={:.6e}".format(
        last['S12'], last['S13'], last['S23']))
    
    print("\n  Volume-averaged STRAIN TENSOR (Voigt):")
    print("    E11={:.6e}  E22={:.6e}  E33={:.6e}".format(
        last['E11'], last['E22'], last['E33']))
    print("    E12={:.6e}  E13={:.6e}  E23={:.6e}".format(
        last['E12'], last['E13'], last['E23']))
    
    print("\n  PRINCIPAL STRESSES:")
    print("    sigma_1={:.6e}  sigma_2={:.6e}  sigma_3={:.6e}".format(
        last['sigma_1'], last['sigma_2'], last['sigma_3']))
    print("    VM={:.6e}  pressure={:.6e}  tresca={:.6e}".format(
        last['sigma_vm'], last['pressure'], last['tresca']))
    
    print("\n  PRINCIPAL STRAINS:")
    print("    eps_1={:.6e}  eps_2={:.6e}  eps_3={:.6e}".format(
        last['epsilon_1'], last['epsilon_2'], last['epsilon_3']))
    
    print("\n  STRAIN ENERGY:  W_density={:.6e}  W_total={:.6e}".format(
        last['W_density'], last['W_total']))
    
    print("\n  ELEMENT STATISTICS:")
    print("    VM: max={:.4e}, mean={:.4e}, std={:.4e}".format(
        last['vm_max'], last['vm_mean'], last['vm_std']))
    
    print("\n  Volume: V_solid={:.6e}, V_RVE={:.6e}".format(
        last['V_solid'], last['V_RVE']))
    
    # ================================================================
    # Write CSV
    # ================================================================
    
    csv_path = os.path.join(Output_Path, 
                             'principal_results_{}.csv'.format(run_id))
    
    with open(csv_path, 'w') as f:
        headers = ['frame', 'time', 
                   'S11', 'S22', 'S33', 'S12', 'S13', 'S23',
                   'E11', 'E22', 'E33', 'E12', 'E13', 'E23',
                   'sigma_1', 'sigma_2', 'sigma_3', 
                   'sigma_vm', 'pressure', 'tresca',
                   'epsilon_1', 'epsilon_2', 'epsilon_3',
                   'W_density', 'W_total',
                   'vm_max', 'vm_mean', 'vm_std',
                   'sp1_max', 'sp1_min', 'sp3_max', 'sp3_min',
                   'press_max', 'press_min']
        f.write(','.join(headers) + '\n')
        
        for r in frame_results:
            vals = [str(r.get(h, 0.0)) for h in headers]
            f.write(','.join(vals) + '\n')
    
    print("\n  Results written to: {}".format(csv_path))
    print("=" * 70 + "\n")
    
    return last


def _applied_strain_from_deck(odb_dir, run_id, suffixes, L):
    """Recover the applied engineering strain as Disp/L from the input deck.

    The decks prescribe a fixed *displacement* on the driving reference point
    (`*Boundary, amplitude=LoadRamp` -> `RP-n, dof, dof, Disp`), so the strain
    the load case actually imposes is Disp/L and therefore depends on the cell
    size. Every load case, tension and shear alike, uses the same magnitude.

    This must not be guessed. A cell of L=0.50 with the standard Disp=0.005
    gives exactly 0.01, so a hardcoded 0.01 is silently right at that one size
    and silently wrong at every other -- it scales the whole tensor by
    (assumed/true) with no error and no warning. Returns None if no deck can be
    read, leaving the caller to decide.
    """
    for sfx in [suffixes[k] for k in sorted(suffixes)]:
        inp = os.path.join(odb_dir, 'Job-{}-{}.inp'.format(run_id, sfx))
        if not os.path.isfile(inp):
            continue
        try:
            armed = False
            for line in open(inp):
                s = line.strip()
                if s.lower().startswith('*boundary') and 'amplitude' in s.lower():
                    armed = True
                    continue
                if armed:
                    if s.startswith('*'):        # section ended without a value
                        armed = False
                        continue
                    parts = [p.strip() for p in s.split(',')]
                    # driving row carries a magnitude: NAME, dof, dof, value
                    if len(parts) >= 4:
                        try:
                            disp = abs(float(parts[3]))
                        except ValueError:
                            continue
                        if disp > 0 and L > 0:
                            return disp / float(L)
        except (IOError, OSError):
            continue
    return None


def extract_elasticity_tensor(odb_dir, output_dir, L, run_id,
                                applied_strain=None,
                                suffixes=None):
    """
    Extract the full 6x6 effective elasticity tensor C_ij from multiple
    load case ODBs.
    
    Requires up to 6 independent load cases:
      - Uniaxial Tension X  (utx): eps_11 = applied_strain
      - Uniaxial Tension Y  (uty): eps_22 = applied_strain
      - Uniaxial Tension Z  (utz): eps_33 = applied_strain
      - Simple Shear S12   (ss12): gamma_12 = applied_strain
      - Simple Shear S13   (ss13): gamma_13 = applied_strain
      - Simple Shear S23   (ss23): gamma_23 = applied_strain
    
    For each load case, the volume-averaged stress tensor gives one
    column of the stiffness matrix:
    
      C_ij = <sigma_i> / eps_j   (no sum on j)
    
    Parameters
    ----------
    odb_dir : str
        Directory containing ODB files.
    output_dir : str
        Directory for output CSV.
    L : float
        RVE side length.
    run_id : str
        RVE identifier (e.g. 'soft_m008').
    applied_strain : float or None
        The engineering strain magnitude applied in each load case. Leave None
        (the default) to read it from the deck as Disp/L, which is the only
        value that is correct at every cell size; pass a float only to override.
    suffixes : dict or None
        Mapping of Voigt column index to job suffix.
        Default: {0:'utx', 1:'uty', 2:'utz', 3:'ss12', 4:'ss13', 5:'ss23'}
    
    Returns
    -------
    C : np.ndarray, shape (6, 6)
        Effective stiffness matrix in Voigt notation.
    S : np.ndarray, shape (6, 6)
        Effective compliance matrix (inverse of C), or zeros if incomplete.
    """
    
    if suffixes is None:
        suffixes = {0: 'utx', 1: 'uty', 2: 'utz',
                    3: 'ss12', 4: 'ss13', 5: 'ss23'}

    # Derive the applied strain from the deck unless explicitly overridden. The
    # whole tensor scales linearly with this number, so getting it wrong scales
    # every modulus by the same factor -- invisible in any ratio, and invisible
    # in the output, which used to print the assumed value as if it were fact.
    strain_src = 'given'
    if applied_strain is None:
        applied_strain = _applied_strain_from_deck(odb_dir, run_id, suffixes, L)
        strain_src = 'from deck (Disp/L)'
        if applied_strain is None:
            raise ValueError(
                "cannot determine the applied strain for {}: no readable "
                "Job-{}-<case>.inp in {}. Pass applied_strain explicitly "
                "(= Disp/L) rather than letting it default -- a wrong value "
                "rescales every modulus silently.".format(run_id, run_id, odb_dir))

    print("\n" + "=" * 70)
    print("ELASTICITY TENSOR EXTRACTION")
    print("=" * 70)
    print("  RVE: {}".format(run_id))
    print("  L: {}".format(L))
    print("  Applied strain: {} [{}]".format(applied_strain, strain_src))
    
    voigt_labels = ['11', '22', '33', '12', '13', '23']
    stress_keys = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
    strain_keys = ['E11', 'E22', 'E33', 'E12', 'E13', 'E23']
    
    C = np.zeros((6, 6))
    stress_tensors = {}
    strain_tensors = {}
    V_solid_stored = None
    V_RVE_stored = None
    available_cases = []
    
    for j, sfx in suffixes.items():
        job_name = 'Job-{}-{}'.format(run_id, sfx)
        odb_path = os.path.join(odb_dir, job_name + '.odb')
        
        if not os.path.isfile(odb_path):
            print("  SKIP column {}: {} not found".format(j, odb_path))
            continue
        
        print("\n  Processing column {} ({}):".format(j, sfx))
        
        result = extract_principals(
            Odb_Path=odb_path,
            Output_Path=output_dir,
            L=L,
            run_id='{}_{}'.format(run_id, sfx),
            last_frame_only=True,
            subtract_mean=False
        )
        
        if result is None:
            print("    ERROR: extraction failed")
            continue
        
        sigma = np.array([result[k] for k in stress_keys])
        epsilon = np.array([result[k] for k in strain_keys])
        
        stress_tensors[sfx] = sigma
        strain_tensors[sfx] = epsilon
        available_cases.append(j)
        
        # Store V_solid and V_RVE from first successful extraction
        if V_solid_stored is None and 'V_solid' in result and 'V_RVE' in result:
            V_solid_stored = result['V_solid']
            V_RVE_stored = result['V_RVE']
        
        # Store stress column (C will be computed as Sigma @ inv(Epsilon) after loop)
        C[:, j] = sigma
        
        print("    <sigma> = [{:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}]".format(*sigma))
        print("    <eps>   = [{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}]".format(*epsilon))
    
    # Compute C = Sigma @ inv(Epsilon)
    # ================================================================
    # HILL-MANDEL CONSISTENT HOMOGENISATION:
    #
    # Stress: sigma_macro = sum(sigma*Ve)/V_RVE (from extract_principals)
    #   -> correctly accounts for void zero-stress contribution
    #
    # Strain: epsilon_macro = applied_strain (from PBC/RP definition)
    #   -> Element-averaged strain misses void deformation and fails
    #      at high porosity. The RP displacement defines the true
    #      macroscopic strain by PBC construction.
    #
    # For each load case j, the applied strain vector has only one
    # non-zero component: epsilon_j = applied_strain (= Disp/L)
    # ================================================================
    E_mat = np.zeros((6, 6))
    
    for j, sfx in suffixes.items():
        if j in available_cases and sfx in strain_tensors:
            # Use applied strain: the PBC enforces exactly this macroscopic strain
            E_mat[j, j] = applied_strain
            # Off-diagonal: read from element field (transverse response is real)
            # These come from extract_principals averaged over V_RVE
            for i in range(6):
                if i != j:
                    E_mat[i, j] = strain_tensors[sfx][i]
    
    if len(available_cases) == 6:
        try:
            C = C @ np.linalg.inv(E_mat)
        except np.linalg.LinAlgError:
            print("\n  WARNING: Strain matrix is singular — falling back to C = sigma/eps_applied")
            for j in available_cases:
                sfx = suffixes[j]
                C[:, j] = stress_tensors[sfx] / applied_strain
    else:
        # Incomplete: fall back to simple division
        for j in available_cases:
            sfx = suffixes[j]
            C[:, j] = stress_tensors[sfx] / applied_strain
    
    # Symmetrise C
    C_sym = 0.5 * (C + C.T)
    asym = np.max(np.abs(C - C_sym))
    
    print("\n" + "=" * 70)
    print("EFFECTIVE STIFFNESS MATRIX C ({} of 6 cases)".format(len(available_cases)))
    print("=" * 70)
    
    print("\n         {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(*voigt_labels))
    for i in range(6):
        row_str = "  {:>2s}  ".format(voigt_labels[i])
        for j in range(6):
            row_str += "{:>12.4e} ".format(C_sym[i, j])
        print(row_str)
    
    print("\n  Max asymmetry: {:.4e}".format(asym))
    
    # Compliance and engineering constants
    S = np.zeros((6, 6))
    if len(available_cases) == 6:
        try:
            S = np.linalg.inv(C_sym)
            
            E1 = 1.0 / S[0, 0]
            E2 = 1.0 / S[1, 1]
            E3 = 1.0 / S[2, 2]
            nu12 = -S[1, 0] * E1
            nu13 = -S[2, 0] * E1
            nu23 = -S[2, 1] * E2
            G12 = 1.0 / S[3, 3]
            G13 = 1.0 / S[4, 4]
            G23 = 1.0 / S[5, 5]
            
            print("\n  ENGINEERING CONSTANTS:")
            print("    E1  = {:.6e}  E2  = {:.6e}  E3  = {:.6e}".format(E1, E2, E3))
            print("    nu12 = {:.4f}        nu13 = {:.4f}        nu23 = {:.4f}".format(nu12, nu13, nu23))
            print("    G12 = {:.6e}  G13 = {:.6e}  G23 = {:.6e}".format(G12, G13, G23))
            
            E_avg = (E1 + E2 + E3) / 3.0
            nu_avg = (nu12 + nu13 + nu23) / 3.0
            G_avg = (G12 + G13 + G23) / 3.0
            G_iso = E_avg / (2 * (1 + nu_avg))
            
            E_aniso = max(abs(E1-E_avg), abs(E2-E_avg), abs(E3-E_avg)) / E_avg * 100
            nu_aniso = max(abs(nu12-nu_avg), abs(nu13-nu_avg), abs(nu23-nu_avg))
            G_aniso = max(abs(G12-G_avg), abs(G13-G_avg), abs(G23-G_avg)) / G_avg * 100
            
            print("\n  ISOTROPY CHECK:")
            print("    E_avg  = {:.6e}  (max deviation: {:.2f}%)".format(E_avg, E_aniso))
            print("    nu_avg = {:.4f}          (max deviation: {:.4f})".format(nu_avg, nu_aniso))
            print("    G_avg  = {:.6e}  (max deviation: {:.2f}%)".format(G_avg, G_aniso))
            print("    G_iso  = E/2(1+nu) = {:.6e}  (vs G_avg: {:.2f}% diff)".format(
                G_iso, abs(G_iso - G_avg) / G_avg * 100))
            
        except np.linalg.LinAlgError:
            print("\n  WARNING: C is singular.")
    else:
        print("\n  NOTE: Only {} of 6 load cases available.".format(len(available_cases)))
        print("  Missing: {}".format([suffixes[j] for j in range(6) if j not in available_cases]))
    
    # Write output
    csv_path = os.path.join(output_dir, 'elasticity_tensor_{}.csv'.format(run_id))
    with open(csv_path, 'w') as f:
        f.write('# Effective stiffness matrix C_ij (Voigt)\n')
        f.write('# RVE: {}, L={}, applied_strain={}\n'.format(run_id, L, applied_strain))
        f.write(',' + ','.join(voigt_labels) + '\n')
        for i in range(6):
            f.write(voigt_labels[i] + ',' + ','.join(
                ['{:.6e}'.format(C_sym[i,j]) for j in range(6)]) + '\n')
        
        f.write('\n# Stress tensors per load case\n')
        f.write('load_case,' + ','.join(['S_'+v for v in voigt_labels]) + '\n')
        for sfx, sigma in stress_tensors.items():
            f.write(sfx + ',' + ','.join(['{:.6e}'.format(s) for s in sigma]) + '\n')
        
        f.write('\n# Strain tensors per load case\n')
        f.write('load_case,' + ','.join(['E_'+v for v in voigt_labels]) + '\n')
        for sfx, epsilon in strain_tensors.items():
            f.write(sfx + ',' + ','.join(['{:.6e}'.format(e) for e in epsilon]) + '\n')
    
    print("\n  Written to: {}".format(csv_path))
    print("=" * 70 + "\n")
    
    return C_sym, S


def extract_principals_batch(odb_dir, output_dir, L, run_ids, 
                              suffixes=['utx', 'ss13', 'ben']):
    """
    Extract principal quantities for multiple RVEs and load cases.
    """
    
    summary = {}
    
    for rid in run_ids:
        for sfx in suffixes:
            job_name = 'Job-{}-{}'.format(rid, sfx)
            odb_path = os.path.join(odb_dir, job_name + '.odb')
            
            if not os.path.isfile(odb_path):
                print("  SKIP: {} not found".format(odb_path))
                continue
            
            label = '{}_{}'.format(rid, sfx)
            print("\n>>> Processing: {}".format(label))
            
            is_bending = sfx.lower() in ['ben', 'bend', 'bending']
            
            try:
                result = extract_principals(
                    Odb_Path=odb_path,
                    Output_Path=output_dir,
                    L=L,
                    run_id=label,
                    last_frame_only=True,
                    subtract_mean=is_bending
                )
                summary[label] = result
            except Exception as e:
                print("  ERROR: {}".format(e))
                summary[label] = None
    
    if summary:
        sum_path = os.path.join(output_dir, 'principal_summary.csv')
        with open(sum_path, 'w') as f:
            headers = ['run_id', 'load_case',
                       'S11', 'S22', 'S33', 'S12', 'S13', 'S23',
                       'E11', 'E22', 'E33', 'E12', 'E13', 'E23',
                       'sigma_1', 'sigma_2', 'sigma_3',
                       'sigma_vm', 'pressure', 'tresca',
                       'epsilon_1', 'epsilon_2', 'epsilon_3',
                       'W_density', 'W_total',
                       'vm_max', 'vm_mean', 'vm_std']
            f.write(','.join(headers) + '\n')
            
            for label, r in summary.items():
                if r is None:
                    continue
                parts = label.rsplit('_', 1)
                rid = parts[0] if len(parts) > 1 else label
                sfx = parts[-1] if len(parts) > 1 else ''
                vals = [rid, sfx] + [str(r.get(h, 0.0)) for h in headers[2:]]
                f.write(','.join(vals) + '\n')
        
        print("\n  Summary written to: {}".format(sum_path))
    
    return summary

# =====================================================================
# STUDY-LEVEL ANALYZERS  (run with plain python3 -- read results CSVs only,
# no Abaqus/ODB). These fold the former standalone analysis scripts into the
# tool so the pipeline needs ONLY SpaX_Standalone / SpaX_GmshPeriodic /
# SpaX_PostProcess. Invoke:  python SpaX_PostProcess.py analyze <name> ...
#   eq19 <results.csv ...>             -- Choi Eq.19 MCST fit + bend/1st-order
#   lengthscale <porous.csv> [base...] -- quadratic l, homogeneous-calibrated
#   homog-calib                        -- homogeneous-cube artifact calibration
#   hybrid <firstorder.csv> <bend.csv> -- small-RVE moduli + big-RVE D_rve -> l
#   rve-study <results.csv> [out_dir]  -- RVE-size convergence (mean/CoV)
# =====================================================================
import math as _math
_R_AVG = 0.04
_D_INC = 2 * _R_AVG

def _a_fnum(s):
    try: return float(s)
    except: return float('nan')

def _a_load(fn):
    return list(csv.DictReader(open(fn, newline=''))) if os.path.isfile(fn) else []

def _a_plt():
    """matplotlib.pyplot (Agg backend) or None if unavailable."""
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None

def _a_converged_moduli():
    """Ensemble E*, nu*, G* from results_lscale.csv (first-order, size-indep)."""
    ls = _a_load('results_lscale.csv')
    E = [_a_fnum(r['E_eff']) for r in ls if _a_fnum(r.get('E_eff', 'nan')) > 0]
    NU = [_a_fnum(r['nu_eff']) for r in ls if _a_fnum(r.get('E_eff', 'nan')) > 0]
    G = [_a_fnum(r['G_eff']) for r in ls if _a_fnum(r.get('E_eff', 'nan')) > 0]
    if not (E and G):
        return float('nan'), 0.348, float('nan')
    return np.mean(E), np.mean(NU), np.mean(G)

def analyze_eq19(csvs):
    """Choi/Lee/Sim Eq.19: E_app(L) = E_classical + 12*mu*(l/L)^2. l from the
    SLOPE of E_app vs 1/L^2 (reference-free). Also the per-RVE bending-vs-first-
    order modulus ratio, separating first-order softening from MCST stiffening."""
    import collections
    Estar, nustar, Gstar = _a_converged_moduli()
    for fn in csvs:
        rows = _a_load(fn)
        by = collections.defaultdict(list)
        for r in rows:
            D = _a_fnum(r.get('D_rve', 'nan')); L = round(_a_fnum(r['L']), 3)
            if not (D > 0): continue
            nu = _a_fnum(r.get('nu_eff', 'nan')); nu = nu if nu == nu else 0.348
            E = _a_fnum(r.get('E_eff', 'nan')); G = _a_fnum(r.get('G_eff', 'nan'))
            Eapp = D * 12.0 / L**4
            Ebm = Eapp * (1.0 - nu**2)
            rat = Ebm / E if E > 0 else float('nan')
            by[L].append((Eapp, Ebm, E, rat, G if G > 0 else Gstar))
        sizes = sorted(by)
        if not sizes:
            print("%s: no valid D_rve rows" % fn); continue
        print("=" * 74); print("FILE: %s" % fn); print("=" * 74)
        print("%-5s %4s %4s | %-9s %-10s %-10s | %s" % (
            "L", "L/d", "n", "E_app(GPa)", "E_bend_mat", "E_1order", "bend/1st"))
        x, y, Gs = [], [], []
        rats = []
        for L in sizes:
            a = np.array(by[L], dtype=float)
            Eapp = np.nanmean(a[:, 0]); Ebm = np.nanmean(a[:, 1])
            Efo = np.nanmean(a[:, 2]); rat = np.nanmean(a[:, 3])
            rats.append(rat); x.append(1.0 / L**2); y.append(Eapp); Gs.append(np.nanmean(a[:, 4]))
            print("%-5.2f %4.1f %4d | %8.3f  %8.3f   %8.3f   | %.3f" % (
                L, L / _D_INC, len(a), Eapp / 1e9, Ebm / 1e9, Efo / 1e9, rat))
        x = np.array(x); y = np.array(y); Gm = np.nanmean(Gs)
        b0, b1 = np.linalg.lstsq(np.vstack([np.ones_like(x), x]).T, y, rcond=None)[0]
        l2 = b1 / (12.0 * Gm)
        print("\nEq.19 fit (mu=G*=%.3g GPa): slope=%+.3e  intercept E0=%.3f GPa" % (
            Gm / 1e9, b1, b0 / 1e9))
        if not _math.isnan(Estar):
            kind = "PLATE" if abs(b0 - Estar / (1 - nustar**2)) < abs(b0 - Estar) else "BEAM"
            print("  intercept vs E*=%.3f / E*plate=%.3f GPa -> RVE bends as a %s" % (
                Estar / 1e9, Estar / (1 - nustar**2) / 1e9, kind))
        if l2 > 0:
            print("  => l = %.4f (l/d=%.2f)  [POSITIVE slope -> MCST stiffening]" % (
                _math.sqrt(l2), _math.sqrt(l2) / _D_INC))
        else:
            print("  => slope<=0 -> l imaginary: NO MCST stiffening")
        print("  bend/1st by size:", "  ".join("%.3f" % r for r in rats),
              "->", "no stiffening beyond first-order softening"
              if max(rats) < 1.03 else "small excess (check size trend: flat=anisotropy, rising=MCST)")
        print()
    return 0

def analyze_lengthscale(porous_fn, baseline_files):
    """Quadratic length scale calibrated by the homogeneous-cube f_quad(L).
    Per RVE l^2=(D_rve - f*Eplate*L^4/12)/(G*L^2), per-RVE moduli when present."""
    import collections
    por = _a_load(porous_fn)
    if not por:
        print("no porous results at %s" % porous_fn); return 1
    Estar, nustar, Gstar = _a_converged_moduli(); Eplate = Estar / (1 - nustar**2)
    if not baseline_files:
        baseline_files = ['results_homog_q_small.csv', 'results_homog_q.csv']
    fq = {}
    for bf in baseline_files:
        for r in _a_load(bf):
            D = _a_fnum(r.get('D_rve', 'nan')); E = _a_fnum(r.get('E_eff', 'nan'))
            nu = _a_fnum(r.get('nu_eff', 'nan')); L = _a_fnum(r['L'])
            if D > 0 and E > 0:
                fq[round(L, 3)] = D / (E / (1 - nu**2) * L**4 / 12.0)
    if not fq:
        print("no homogeneous baseline (%s)" % ", ".join(baseline_files)); return 1
    print("Quadratic homogeneous baseline f_quad(L):")
    for L in sorted(fq): print("  L=%.2f (L/d=%.1f) f=%.4f" % (L, L / _D_INC, fq[L]))
    bysize = collections.defaultdict(list); nper = 0
    for r in por:
        D = _a_fnum(r.get('D_rve', 'nan')); L = round(_a_fnum(r['L']), 3)
        if not (D > 0): continue
        f = fq.get(L) or fq[min(fq, key=lambda k: abs(k - L))]
        Er = _a_fnum(r.get('E_eff', 'nan')); nur = _a_fnum(r.get('nu_eff', 'nan'))
        Gr = _a_fnum(r.get('G_eff', 'nan'))
        if Er > 0 and Gr > 0 and nur == nur:
            Ep, Gv = Er / (1 - nur**2), Gr; nper += 1
        else:
            Ep, Gv = Eplate, Gstar
        Dcl = f * Ep * L**4 / 12.0
        l2 = (D - Dcl) / (Gv * L**2); l = _math.sqrt(l2) if l2 > 0 else 0.0
        bysize[L].append((l, D / Dcl))
    print("\nPer-RVE moduli used for %d/%d RVEs.\n" % (
        nper, sum(len(v) for v in bysize.values())))
    sizes = sorted(bysize); Sm1 = []; inv = []
    print("%-6s %5s %4s | l (mean std CoV%%) | S-1" % ("L", "L/d", "n"))
    for L in sizes:
        a = np.array(bysize[L]); lm = a[:, 0].mean()
        lsd = a[:, 0].std(ddof=1) if len(a) > 1 else 0.0
        s1 = a[:, 1].mean() - 1; Sm1.append(s1); inv.append(1.0 / L**2)
        print("%-6.2f %5.1f %4d | %.4f %.4f %5.1f | %+.4f" % (
            L, L / _D_INC, len(a), lm, lsd, lsd / lm * 100 if lm > 0 else float('nan'), s1))
    if len(sizes) >= 2:
        x = np.array(inv); yv = np.array(Sm1)
        slope = np.dot(x, yv) / np.dot(x, x)
        l2 = slope * 12 * Gstar / Eplate
        print("\n(S-1)~slope/L^2: slope=%.3e -> l_fit=%.4f%s" % (
            slope, _math.sqrt(l2) if l2 > 0 else 0.0,
            "" if l2 > 0 else "  (slope<=0: no length scale)"))
    return 0

def analyze_homog_calib():
    """Homogeneous-cube calibration: f = D_rve/(E_plate L^4/12) for true-l=0
    cubes -- the cube-vs-thin-plate + (linear) locking artifact, by size & mesh."""
    def dr(r):
        D = _a_fnum(r.get('D_rve', 'nan')); E = _a_fnum(r.get('E_eff', 'nan'))
        nu = _a_fnum(r.get('nu_eff', 'nan')); L = _a_fnum(r['L'])
        return D / (E / (1 - nu**2) * L**4 / 12.0) if (D > 0 and E > 0) else None
    print("Homogeneous-cube artifact factor f (true l=0):")
    for fn in ('results_homog_small.csv', 'results_homog.csv',
               'results_homog_q_small.csv', 'results_homog_q.csv'):
        rows = _a_load(fn)
        if not rows: continue
        print("--- %s ---" % fn)
        for r in sorted(rows, key=lambda x: _a_fnum(x['L'])):
            f = dr(r)
            print("  L=%.2f (L/d=%.1f)  f=%s" % (
                _a_fnum(r['L']), _a_fnum(r['L']) / _D_INC,
                "%.4f" % f if f else "UNSOLVED"))
    mc = _a_load('results_homog_meshconv.csv')
    if mc:
        print("--- mesh convergence (results_homog_meshconv.csv) ---")
        for r in sorted(mc, key=lambda x: -_a_fnum(x['L_mesh'])):
            f = dr(r)
            print("  L_mesh=%s f=%s" % (r['L_mesh'], "%.4f" % f if f else "UNSOLVED"))
    return 0

def analyze_hybrid(fo_fn, bend_fn, out='.'):
    """Hybrid l: converged first-order moduli (small RVEs) + per-big-RVE D_rve."""
    import collections
    FO_MIN_LD = 5.0
    fo = _a_load(fo_fn); bd = _a_load(bend_fn)
    E = [_a_fnum(r['E_eff']) for r in fo
         if _a_fnum(r['L']) / _D_INC >= FO_MIN_LD and _a_fnum(r.get('E_eff', 'nan')) > 0]
    NU = [_a_fnum(r['nu_eff']) for r in fo
          if _a_fnum(r['L']) / _D_INC >= FO_MIN_LD and _a_fnum(r.get('E_eff', 'nan')) > 0]
    G = [_a_fnum(r['G_eff']) for r in fo
         if _a_fnum(r['L']) / _D_INC >= FO_MIN_LD and _a_fnum(r.get('G_eff', 'nan')) > 0]
    if not (E and NU and G):
        print("no first-order E/nu/G at L/d>=%.1f" % FO_MIN_LD); return 1
    Estar, nustar, Gstar = np.mean(E), np.mean(NU), np.mean(G)
    print("Converged first-order: E*=%.4g GPa nu*=%.3f G*=%.4g GPa\n" % (
        Estar / 1e9, nustar, Gstar / 1e9))
    bysize = collections.defaultdict(list)
    for r in bd:
        D = _a_fnum(r.get('D_rve', 'nan')); L = _a_fnum(r['L'])
        if not (D > 0): continue
        Dcl = Estar / (1 - nustar**2) * L**4 / 12.0
        l2 = (D - Dcl) / (Gstar * L**2); l = _math.sqrt(l2) if l2 > 0 else 0.0
        bysize[round(L, 4)].append((l, D / Dcl))
    print("%-6s %5s %4s | l (mean std CoV%%) | D_ratio" % ("L", "L/d", "n"))
    for L in sorted(bysize):
        a = np.array(bysize[L]); m = a[:, 0].mean()
        s = a[:, 0].std(ddof=1) if len(a) > 1 else 0.0
        print("%-6.2f %5.1f %4d | %.4f %.4f %5.1f | %.3f" % (
            L, L / _D_INC, len(a), m, s, s / m * 100 if m > 0 else float('nan'),
            a[:, 1].mean()))
    return 0

def analyze_rve_study(res_fn, out='.'):
    """RVE-size convergence: mean and seed CoV of E_eff/G_eff/l/D_ratio vs L/d."""
    COV_TOL, PLATEAU_TOL = 0.02, 0.02
    METRICS = [("E_eff", "E_eff", 1e9), ("G_eff", "G_eff", 1e9),
               ("l", "l", 1.0), ("D_ratio", "D_rve/D_class", 1.0)]
    rows = _a_load(res_fn)
    bysize = {}
    for r in rows:
        bysize.setdefault(round(_a_fnum(r['L']), 4), []).append(r)
    sizes = sorted(bysize)
    if not sizes:
        print("no rows in %s" % res_fn); return 1
    stats = {m[0]: {} for m in METRICS}
    print("%-6s %5s %4s | %s" % ("L", "L/d", "n",
          " | ".join("%-20s" % m[1] for m in METRICS)))
    for L in sizes:
        grp = bysize[L]; line = "%-6.2f %5.1f %4d | " % (L, L / _D_INC, len(grp))
        for key, lbl, sc in METRICS:
            vals = np.array([_a_fnum(r.get(key, 'nan')) for r in grp])
            vals = vals[~np.isnan(vals)] / sc
            if len(vals) == 0:
                stats[key][L] = (float('nan'),) * 3; line += "%20s | " % "-"; continue
            mean = vals.mean(); std = vals.std(ddof=1) if len(vals) > 1 else 0.0
            cov = std / mean * 100 if mean != 0 else float('nan')
            stats[key][L] = (mean, std, cov)
            line += "%7.4g %6.3g %4.1f%% | " % (mean, std, cov)
        print(line)
    Lmax = sizes[-1]
    print("\nRepresentative size (CoV<%.0f%% AND mean within %.0f%% of largest):" % (
        COV_TOL * 100, PLATEAU_TOL * 100))
    for key, lbl, sc in METRICS:
        ref = stats[key][Lmax][0]; rec = None
        for L in sizes:
            m, s, cov = stats[key][L]
            if _math.isnan(m) or ref == 0: continue
            if cov / 100 <= COV_TOL and abs(m - ref) / abs(ref) <= PLATEAU_TOL:
                rec = L; break
        print("  %-14s : %s" % (lbl, "L/d>=%.1f" % (rec / _D_INC) if rec
                                else "not converged in range"))
    return 0

def _analyze(args):
    if not args:
        print("analyze subcommands: eq19 | lengthscale | homog-calib | hybrid | rve-study")
        return 1
    name, rest = args[0], args[1:]
    if name == 'eq19':
        return analyze_eq19(rest or ['results.csv'])
    if name == 'lengthscale':
        if not rest: print("usage: analyze lengthscale <porous.csv> [baseline.csv ...]"); return 1
        return analyze_lengthscale(rest[0], rest[1:])
    if name == 'homog-calib':
        return analyze_homog_calib()
    if name == 'hybrid':
        if len(rest) < 2: print("usage: analyze hybrid <firstorder.csv> <bending.csv>"); return 1
        return analyze_hybrid(rest[0], rest[1])
    if name == 'rve-study':
        if not rest: print("usage: analyze rve-study <results.csv> [out_dir]"); return 1
        return analyze_rve_study(rest[0], rest[1] if len(rest) > 1 else '.')
    print("unknown analyze subcommand: %s" % name); return 1


if __name__ == '__main__':
    # Study-level analyzers (plain python3, read CSVs only -- no Abaqus).
    if len(sys.argv) >= 2 and sys.argv[1] == 'analyze':
        sys.exit(_analyze(sys.argv[2:]))

    # Principal / full-tensor / elasticity-tensor extraction (abaqus python).
    if len(sys.argv) >= 2 and sys.argv[1] == 'principal':
        if len(sys.argv) < 4:
            print("Usage: abaqus python SpaX_PostProcess.py principal <odb> <out> [L] [run_id]")
            sys.exit(1)
        extract_principals(sys.argv[2], sys.argv[3],
                           float(sys.argv[4]) if len(sys.argv) > 4 else 0.0,
                           sys.argv[5] if len(sys.argv) > 5 else 'RVE')
        sys.exit(0)
    if len(sys.argv) >= 2 and sys.argv[1] == 'elasticity':
        if len(sys.argv) < 6:
            print("Usage: abaqus python SpaX_PostProcess.py elasticity "
                  "<odb_dir> <out> <L> <run_id> [applied_strain]")
            print("       applied_strain is read from the deck as Disp/L if omitted;")
            print("       pass it only to override.")
            sys.exit(1)
        _eps = float(sys.argv[6]) if len(sys.argv) > 6 else None
        extract_elasticity_tensor(sys.argv[2], sys.argv[3], float(sys.argv[4]),
                                  sys.argv[5], applied_strain=_eps)
        sys.exit(0)

    # Merge mode (no Abaqus needed): combine per-RVE partial CSVs.
    if len(sys.argv) >= 2 and sys.argv[1] == '--merge':
        if len(sys.argv) < 4:
            print("Usage: python SpaX_PostProcess.py --merge <parts_dir> <output_csv>")
            sys.exit(1)
        merge_partials(sys.argv[2], sys.argv[3])
        sys.exit(0)

    if len(sys.argv) < 3:
        print("Usage: abaqus python SpaX_PostProcess.py <csv_path> <odb_dir> [output_csv] [index]")
        print("       (index = 1-based RVE row; omit to process all RVEs in one pass)")
        print("       python SpaX_PostProcess.py --merge <parts_dir> <output_csv>")
        print("       python SpaX_PostProcess.py analyze <eq19|lengthscale|homog-calib|hybrid|rve-study> ...")
        sys.exit(1)

    csv_path = sys.argv[1]
    odb_dir = sys.argv[2]
    output_csv = sys.argv[3] if len(sys.argv) > 3 else 'results.csv'
    only_index = int(sys.argv[4]) if len(sys.argv) > 4 else None

    results = run(csv_path, odb_dir, output_csv, only_index=only_index)
    if only_index is None:
        print_summary_table(results)
