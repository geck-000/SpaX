# -*- coding: utf-8 -*-
"""
Spatium_PrincipalExtract.py — Extract principal stresses, strains, 
                               strain energy, full tensors, and 
                               elasticity tensor from RVE ODBs
===================================================================

Produces volume-averaged and per-element principal quantities for
multi-axial characterisation of RVE response.

Usage from Abaqus CAE console:
    import Spatium_PrincipalExtract as PE
    
    # Single ODB extraction
    PE.extract_principals(
        Odb_Path=r'C:\\SIMULIA\\temp\\Job-myRVE-utx.odb',
        Output_Path=r'C:\\SIMULIA\\temp',
        L=0.2,
        run_id='soft_m008'
    )
    
    # Full elasticity tensor from multiple load cases
    PE.extract_elasticity_tensor(
        odb_dir=r'C:\\SIMULIA\\temp',
        output_dir=r'C:\\SIMULIA\\temp',
        L=0.2,
        run_id='soft_m008',
        applied_strain=0.01
    )
"""

from __future__ import print_function
import os
import math
import numpy as np

from odbAccess import openOdb


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


def _get_data(field_value):
    """Get data from a field value, handling both single and double precision."""
    try:
        return field_value.data
    except:
        return field_value.dataDouble


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
    
    odb = openOdb(path=Odb_Path, readOnly=True)
    
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


def extract_elasticity_tensor(odb_dir, output_dir, L, run_id,
                                applied_strain=0.01,
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
    applied_strain : float
        The engineering strain magnitude applied in each load case.
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
    
    print("\n" + "=" * 70)
    print("ELASTICITY TENSOR EXTRACTION")
    print("=" * 70)
    print("  RVE: {}".format(run_id))
    print("  Applied strain: {}".format(applied_strain))
    
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
