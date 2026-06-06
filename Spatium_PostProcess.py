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
from odbAccess import openOdb


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


def extract_first_order(odb_path, s_comp, eng_strain, L):
    """
    Extract E or G and nu from an ODB using direct odbAccess.
    """
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
            s_sub = stress_field.getSubset(region=inst)
            v_sub = volume_field.getSubset(region=inst)
            vol_dict = {v.elementLabel: _get_data(v) for v in v_sub.values}
            
            # Strain field
            strain_dict = {}
            if strain_key:
                e_sub = frame.fieldOutputs[strain_key].getSubset(region=inst)
                for ev in e_sub.values:
                    strain_dict[ev.elementLabel] = _get_data(ev)
            
            for sv in s_sub.values:
                el = sv.elementLabel
                if el not in vol_dict:
                    continue
                vol = vol_dict[el]
                s_data = _get_data(sv)
                total_stress_vol += s_data[s_idx] * vol
                total_vol += vol
                
                # Axial strain (element-averaged)
                if el in strain_dict:
                    e_data = strain_dict[el]
                    total_strain_axial_vol += e_data[s_idx] * vol
                    
                    # Transverse strains for Poisson's ratio
                    if not is_shear:
                        if s_comp == 'S11':
                            total_strain_trans1_vol += e_data[1] * vol
                            total_strain_trans2_vol += e_data[2] * vol
                        elif s_comp == 'S22':
                            total_strain_trans1_vol += e_data[0] * vol
                            total_strain_trans2_vol += e_data[2] * vol
                        elif s_comp == 'S33':
                            total_strain_trans1_vol += e_data[0] * vol
                            total_strain_trans2_vol += e_data[1] * vol
        
        # Hill-Mandel consistent:
        # Stress: sum(sigma*Ve) / V_RVE
        # Axial strain: RP-based (frameValue * eng_strain)
        # Transverse strain: element-averaged / V_solid (for Poisson only)
        sigma_macro = total_stress_vol / V_RVE
        eps_macro = frame.frameValue * eng_strain
        
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
    
    # Process last frame
    last_frame = step.frames[-1]
    stress_field = last_frame.fieldOutputs['S']
    volume_field = last_frame.fieldOutputs['EVOL']
    
    sigma_vol = np.zeros(6)
    moment_vol = 0.0
    total_vol = 0.0
    
    for inst_name, inst in solid_instances:
        s_sub = stress_field.getSubset(region=inst)
        v_sub = volume_field.getSubset(region=inst)
        vol_dict = {v.elementLabel: _get_data(v) for v in v_sub.values}
        
        for sv in s_sub.values:
            el = sv.elementLabel
            key = (inst_name, el)
            if el not in vol_dict or key not in elem_centroid:
                continue
            vol = vol_dict[el]
            s_data = _get_data(sv)
            centroid = elem_centroid[key]
            z_rel = centroid[coord_idx] - coord_bar
            
            for c in range(min(6, len(s_data))):
                sigma_vol[c] += s_data[c] * vol
            moment_vol += s_data[s_idx] * z_rel * vol
            total_vol += vol
    
    kappa_actual = Kappa * last_frame.frameValue
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


def run(csv_path, odb_dir, output_csv='postprocess_results.csv'):
    """Main post-processing pipeline."""
    
    params_list = read_csv(csv_path)
    print("\n" + "=" * 70)
    print("SPATIUM POST-PROCESSING PIPELINE")
    print("=" * 70)
    print("  CSV: {}".format(csv_path))
    print("  ODB dir: {}".format(odb_dir))
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


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: abaqus python Spatium_PostProcess.py <csv_path> <odb_dir> [output_csv]")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    odb_dir = sys.argv[2]
    output_csv = sys.argv[3] if len(sys.argv) > 3 else 'results.csv'
    
    results = run(csv_path, odb_dir, output_csv)
    print_summary_table(results)
