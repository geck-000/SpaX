# hpc/ — HPC batch scripts (CSC Roihu / Slurm)

Slurm scripts used to run the campaign on CSC Roihu. They are archival examples
of the generate → solve → post-process pipeline at cluster scale; paths and
account/project fields are specific to that environment and should be adapted
before reuse.

## Site-specific fields

Two placeholders stand in for the allocation these were run under. Neither is a
working value — set them for your own site:

| Placeholder | Meaning | How to set |
|---|---|---|
| `project_XXXXXX` | Slurm account / allocation | edit the `#SBATCH --account=` line, or override per submission with `sbatch --account=<acct> …`, which takes precedence over the directive |
| `/scratch/project_XXXXXX/test_rve` | shared working directory on the cluster | export `WORKDIR=…` before running; every script now honours it |

`PYTHONUSERBASE` is likewise overridable where the scripts set it. The scripts
otherwise assume the core modules (`SpaX_Standalone.py`, `SpaX_PostProcess.py`)
and the decks are staged in `WORKDIR`.

| File | Role |
|------|------|
| `submit_firstorder.sh` | Slurm array: solve the first-order (uniaxial/shear) decks. |
| `submit_coltensor.sh` | Solve the full-column 6×6 tensor decks. |
| `submit_brineK.sh` | Solve the brine `K(T)` study [#5] decks. |
| `submit_nlgeom.sh` | Solve the large-deformation study [#8] decks. |
| `submit_tilt.sh` | Solve the channel-tilt study [#6] decks. |
| `submit_colseeds.sh` | Solve the seeded C-shape column (`rve_colseeds.csv`, 100 decks) for the depth-profile scatter envelopes, then postprocess to `results_colseeds.csv`. |
| `rerun_fieldseeds_colseeds.sh` | Staged controller for the ten-packing re-run after the seeding fix: solves the 300 fieldseeds decks in two waves, then generates and solves the 100 colseeds-extra decks. Re-submits itself per `$STAGE`, so it must be staged in `WORKDIR`. Chained because the `small` partition caps submitted jobs at 200. |
| `postprocess_firstorder.sh` | `abaqus python SpaX_PostProcess.py` over the solved first-order ODBs. |
| `postprocess_coltensor.sh` | Extract the per-slice 6×6 elasticity tensors. |
| `postprocess_nlgeom.sh` | Extract nlgeom reference-point curves. |

### Earlier campaign scripts

Restored from the repository root, where they had been left untracked:

| File | Role |
|------|------|
| `csc_solve_array.sh` | The general Slurm solve array used for most campaigns. |
| `bend_01_generate.sh` … `bend_04_merge.sh`, `bend_submit.sh` | The four-stage bending (second-order) pipeline and its submitter. |
| `run_channels_q.sh`, `run_porous_q_full.sh` | Quadratic-element channel and porous sweeps. |
| `submit_failure.sh`, `postprocess_failure.sh` | Failure-onset campaign. |
| `submit_si2nd_l400.sh`, `postprocess_si2nd_l400.sh` | The `L=0.40` second-order set, the heavy part of the 3-size MCST null. |

> **Trap.** `csc_solve_array.sh` skips a deck when its `.odb` merely *exists*,
> not when the solve completed. A job killed at walltime leaves a truncated ODB,
> so a naive resubmit silently skips exactly the jobs that failed. Identify
> failures by the absence of an `Abaqus exit:` line in the log and delete those
> ODBs before resubmitting. See `../docs/RUNBOOK.md` §5.

**Note.** After the CSC v2026_03 software-stack change, `module load abaqus`
fails in batch; these scripts source a snapshot env (`~/abaqus_env.sh`) instead.
They assume the core modules (`SpaX_Standalone.py`, `SpaX_PostProcess.py`) and
the decks are staged in the job working directory.

## CSC Roihu access

The scripts submit from a Roihu login node. CSC authenticates with short-lived
SSH **certificates**: a login tool generates a keypair and has CSC's CA sign the
public half, so `ssh` needs *both* the private key and the signed certificate.

Define a `roihu` host in `~/.ssh/config` pointing at `~/.ssh/csc_id` +
`~/.ssh/csc_id-cert.pub`, with your CSC username as `User`. Two things are not
in the repository and must be present locally before anything can be submitted:

1. **The private key.** A certificate on its own cannot authenticate. Re-running
   the CSC login tool writes `csc_id` and `csc_id-cert.pub` together; a
   certificate downloaded by itself from the web UI is not enough. Certificates
   expire after 24 h, so this is a routine per-session step.
2. **Network access.** CSC restricts SSH to allowlisted networks. From outside
   them `roihu.csc.fi:22` times out (while, say, `github.com:22` connects
   normally, which is the quick way to tell a CSC allowlist problem from a local
   firewall). Connect through the institutional VPN first.

Check both with:

```bash
ssh-keygen -L -f ~/.ssh/csc_id-cert.pub   # principal + validity window
ssh -o BatchMode=yes roihu true           # should not time out
```
