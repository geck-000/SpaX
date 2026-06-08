"""
Spatium_PostProcess.py
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

    Returns (labels: int ndarray (n,), data: float ndarray (n, ncomp)).
    Assumes one value per element (1 integration point, as for C3D4/C3D4H);
    that matches the volume-weighting the callers already do with whole-element
    EVOL.
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
            return np.concatenate(labs), np.concatenate(dats, axis=0)

    # Fallback: per-element .values (slow path).
    labs, dats = [], []
    for v in sub.values:
        labs.append(v.elementLabel)
        dats.append(np.asarray(_get_data(v), dtype=float).reshape(-1))
    if not labs:
        return np.array([], dtype=int), np.zeros((0, 1))
    return np.asarray(labs, dtype=int), np.vstack(dats)


def _aligned(values_map, labels):
    """Vectorized lookup: array of values_map[l] for l in labels, 0.0 if absent."""
    return np.fromiter((values_map.get(l, 0.0) for l in labels),
                       dtype=float, count=len(labels))


# Maps the driven stress component to the reference-point node set and the
# 0-based DOF on it that carries the imposed macro deformation. Normal modes
# drive RP-1/2/3 on their normal DOF; shear modes drive RP-4 on the shear DOF
# (ss12->2, ss13->1, ss23->2, matching Spatium_Standalone.shear_config).
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
        })
    
    odb.close()
    
    results = {}
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
    print("SPATIUM POST-PROCESSING PIPELINE")
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
        
        # ---- First-order: UT ----
        odb1 = os.path.join(odb_dir, 'Job-{}-{}.odb'.format(run_id, mode_short(Mode)))
        if os.path.isfile(odb1):
            print("  Extracting {} ({})...".format(mode_short(Mode), stress_component(Mode)))
            try:
                eng_strain = Disp / L
                r1 = extract_first_order(odb1, stress_component(Mode), eng_strain, L)
                row['E_eff'] = r1.get('E_eff', '')
                row['nu_eff'] = r1.get('nu_eff', '')
                print("    E = {:.4e}, nu = {}".format(
                    r1.get('E_eff', 0), r1.get('nu_eff', '')))
            except Exception as e:
                print("    ERROR: {}".format(e))
                row['E_eff'] = 'ERROR'
        else:
            print("  [SKIP] {} not found".format(os.path.basename(odb1)))
            row['E_eff'] = 'MISSING'
        
        # ---- First-order: Shear ----
        if Mode2:
            odb2 = os.path.join(odb_dir, 'Job-{}-{}.odb'.format(run_id, mode_short(Mode2)))
            if os.path.isfile(odb2):
                print("  Extracting {} ({})...".format(mode_short(Mode2), stress_component(Mode2)))
                try:
                    eng_strain2 = Disp2 / L
                    r2 = extract_first_order(odb2, stress_component(Mode2), eng_strain2, L)
                    row['G_eff'] = r2.get('G_eff', r2.get('E_eff', ''))
                    print("    G = {:.4e}".format(float(row['G_eff'])))
                except Exception as e:
                    print("    ERROR: {}".format(e))
                    row['G_eff'] = 'ERROR'
            else:
                print("  [SKIP] {} not found".format(os.path.basename(odb2)))
                row['G_eff'] = 'MISSING'
        
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


if __name__ == '__main__':
    # Merge mode (no Abaqus needed): combine per-RVE partial CSVs.
    if len(sys.argv) >= 2 and sys.argv[1] == '--merge':
        if len(sys.argv) < 4:
            print("Usage: python Spatium_PostProcess.py --merge <parts_dir> <output_csv>")
            sys.exit(1)
        merge_partials(sys.argv[2], sys.argv[3])
        sys.exit(0)

    if len(sys.argv) < 3:
        print("Usage: abaqus python Spatium_PostProcess.py <csv_path> <odb_dir> [output_csv] [index]")
        print("       (index = 1-based RVE row; omit to process all RVEs in one pass)")
        print("       python Spatium_PostProcess.py --merge <parts_dir> <output_csv>")
        sys.exit(1)

    csv_path = sys.argv[1]
    odb_dir = sys.argv[2]
    output_csv = sys.argv[3] if len(sys.argv) > 3 else 'results.csv'
    only_index = int(sys.argv[4]) if len(sys.argv) > 4 else None

    results = run(csv_path, odb_dir, output_csv, only_index=only_index)
    if only_index is None:
        print_summary_table(results)
