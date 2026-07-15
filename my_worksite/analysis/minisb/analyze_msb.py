#!/usr/bin/env python3
"""
analyze_msb.py
======================
Scientific analysis script for the Mini Santa Barbara test (Nyz).

Produces:
  - Thermal history across redshift for maximum and density weighted temperature
  - 

Usage:
    python analyze_msb.py --plotfiles /path/to/output/plt* --save /path/to/save
"""

import argparse
import glob
import os
import sys
import yaml

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from natsort import natsorted
from unyt import cm, km, s, Mpc
from datetime import datetime

import yt
yt.set_log_level(40)  # suppress yt output except critical errors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Analyze Mini Santa Barbara Nyx run")
    p.add_argument("--plotfiles", default=None,
                   help="Glob pattern or superior directory for plotfiles, logs and headers (default: auto-detect)")
    p.add_argument("--save", default="./outputs",
                   help="Directory to save plots and animations")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------
def create_runlog_df(runlog_path):
    """ Return the pandas dataframe of the target columns of the runlog """

    # List inizialization from runlog headers
    target_cols =["nstep", "time", "dt","z", "a", "T_max", "T_at_rho", "T_rho"]   
    headers = {head: [] for head in target_cols}

    if not os.path.exists(runlog_path):
        print(f"Error: {runlog_path} not found")
        return 1
    
    with open(runlog_path, 'r') as f:
        for i, line in enumerate(f):    
            line = line.strip()                           # clean tabulations and spaces
            line = line.replace("T @ <rho>", "T_at_rho")  # usable format
            line = line.replace("<T>_rho", "T_rho")  
            line = line.replace("<T>_V", "T_V")  
            line = line.replace("T(21cm)", "T21")  

            if not line:
                continue

            columns = line.split() # list of the row, each entry is a column
            if i == 0:   # header
                columns.pop(0)     # remove the first # in the header
                idx_map = {name: columns.index(name) for name in columns}  # dic: name = header, value = idx pos in runlog
            else:        # data
                for name in target_cols:
                    idx = idx_map[name]     # get the index of the header
                    value = float(columns[idx])
                    headers[name].append(value)

    return pd.DataFrame(headers) # convert into pandas dataframe
    
            
# ---------------------------------------------------------------------------
# Output directory setup
# ---------------------------------------------------------------------------
def make_dirs(base):
    subdirs = ["ThermalEvolution", "Gas", "Gravity"]
    for d in subdirs:
        os.makedirs(os.path.join(base, d), exist_ok=True)


# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    save_path = args.save
    make_dirs(save_path)

    if args.plotfiles:
        plotfiles = natsorted(glob.glob(os.path.join(args.plotfiles, "plt*")))                                       
    else:  # try to auto-detect in common output locations
        candidates = glob.glob("/data/mfulghieri/Nyx/outputs/mini_sb/msb_*")
        if not candidates:
            print("No plotfiles found. Use --plotfiles."); sys.exit(1)
        run_dir = natsorted(candidates)[-1]   # most recent run
        plotfiles = natsorted(glob.glob(os.path.join(run_dir, "plt*")))
        print(f"Auto-detected run: {run_dir}")
        
    plotfiles = [p for p in plotfiles if not p.endswith(".old")]

    if not plotfiles:
        print("No plotfiles found!"); sys.exit(1)
    print(f"Found {len(plotfiles)} plotfiles.")

    # ---- Runlog statistical analysis ----
    runlog_path = os.path.join(args.plotfiles, "runlog")
    df = create_runlog_df(runlog_path)   # dataframe of selected runlog columns

    print("\n" + "="*60)
    print("      SANTA BARBARA TEST MINI - RUNLOG DIAGNOSTICS      ")
    print("="*60)

    # Simulation limits
    row_start = df.iloc[0]
    row_end = df.iloc[-1]
    print(f"     COSMOLOGICAL LIMITS:")
    print(f"   - Initial redshift (z): {row_start['z']:.2f}  ->  Final: {row_end['z']:.2f}")
    print(f"   - Scale factor (a)    :  {row_start['a']:.4f}  ->  Final: {row_end['a']:.4f}")
    print(f"   - Total steps (nstep) :  {int(row_end['nstep'])}")

    # Thermal peak
    idx_shock = df['T_max'].idxmax()
    row_shock = df.loc[idx_shock]
    print(f"\n  MAXIMUM COLLAPSE MOMENT (SHOCK):")
    print(f"   - Step        :  {int(row_shock['nstep'])}")
    print(f"   - Redshift (z):  {row_shock['z']:.3f}")
    print(f"   - T_max       :  {row_shock['T_max']:.3e} K")
    print(f"   - T_rho       :  {row_shock['T_rho']:.3e} K")

    # Timestep
    dt_min = df['dt'].min()
    dt_max = df['dt'].max()
    print(f"\n   TIMESTEP (dt):")
    print(f"   - Minimum dt (max effort):  {dt_min:.3e}")
    print(f"   - Maximum dt (min effort):  {dt_max:.3e}")

    # Summary
    print("\n  THERMODYNAMIC SUMMARY:")
    thermal_terms = df[['T_max', 'T_at_rho', 'T_rho']]
    print(thermal_terms.describe().loc[['min', 'mean', 'max']].to_string(formatters={  # .tostring() to fancy formatting
        'T_max': '{:,.2e}'.format,
        'T_at_rho': '{:,.2e}'.format,
        'T_rho': '{:,.2e}'.format
    }))
    print("="*60 + "\n")

    # ---- Thermal Evolution Plot ----
    fig_t, ax_t = plt.subplots(figsize=(10, 5))
    ax_t.plot(df['z'], df['T_max'], lw=2.5, color="#1f77b4", label=r"$T_{\max}$ (Local Peak)")
    ax_t.plot(df['z'], df['T_rho'], lw=2, ls="--", color="crimson", label=r"$\langle T \rangle_{\rho}$ (Density Weighted)")
    ax_t.set_title("Santa Barbara Test Mini | Thermal Evolution & Shock History", fontweight="bold", fontsize=12)
    ax_t.set_xlabel("Redshift ($z$)", fontsize=11)
    ax_t.set_ylabel("Temperature ($K$)", fontsize=11)
    ax_t.set_yscale('log')
    ax_t.set_xlim(df['z'].max(), df['z'].min()) # Iinvert z x-ax
    ax_t.set_ylim(df['T_rho'].min() * 0.5, df['T_max'].max() * 3.0)
    ax_t.legend(loc="upper left")
    ax_t.grid(True, ls=":", alpha=0.5)

    fig_t.tight_layout()
    fig_t.savefig(os.path.join(save_path, "ThermalEvolution", "runlog_thermal_history.png"), dpi=200)
    plt.close(fig_t)



    # --- Storage for history ---
    times    = []
    z_values = []
    a_values = []

    for i, plt_path in enumerate(plotfiles):
        ds  = yt.load(plt_path, hint='Nyx')   # both Nyx and amrex format possible
        ad  = ds.all_data()

        # ---- Metadata ----
        time = ds.current_time
        z_now = ds.current_redshift   # already given in the runlog
        a_now = 1 / (1 + z_now)
        times.append(time)
        z_values.append(z_now)
        a_values.append(a_now)

        # Physical box size
        Lx = float(ds.domain_width[0].to("Mpc").v)
        Ly = float(ds.domain_width[1].to("Mpc").v)


        # ---- 2D density slice ----
        fields = [('boxlib', 'density'), ('boxlib', 'rho_e'), ('boxlib', 'Temp'), 
                  ('boxlib', 'xmom'), ('boxlib', 'ymom'), ('boxlib', 'zmom'), 
                  ('boxlib', 'magvel'), ('boxlib', 'phi_grav'), ('boxlib', 'particle_mass_density')]
        # ('boxlib', 'density'): gas mass density. ('boxlib', 'rho_e'): gas internal energy density.
        # ('boxlib', 'Temp'): Gas temperature. ('boxlib', 'xmom'/'ymom'/'zmom'): Gas momentum.
        # ('boxlib', 'magvel'): Magnitude of gas velocity (|v|). 
        # ('boxlib', 'phi_grav'): Total gravitational potential (phi). ('boxlib', 'particle_mass_density'):
        # DM mass density projected onto the grid. 
        slc = yt.SlicePlot(ds, 'z', fields, center='c')   
        
        # ---- Density fixed resolution buffer (FRB) matrix (res x res) ----
        res = 1024
        slc.set_buff_size(res)                                # pixel resolution of the image
        frb = slc.frb[('boxlib', 'density')].v                # extract the numerical matrix (res, res)
        gas_density_log = np.log10(np.clip(frb, 1e-32, None)) # log scale, lower clip to avoid log <= 0, no upper clip
      
        fig_gdens_frb, ax_gdens_frb = plt.subplots(figsize=(8, 8), facecolor='black')
        ax_gdens_frb.set_facecolor('black')   # black plot background (visible if missing data)
        im = ax_gdens_frb.imshow(gas_density_log.T, origin="lower", 
                           extent=[0, Lx, 0, Ly],
                           cmap="magma",             
                           interpolation="bilinear") # gradient effect between pixels
        ax_gdens_frb.set_title(f"Cosmic Gas Density  |  $z = {z_now:.2f}$", color="white", fontsize=13, fontweight='bold', pad=15)
        ax_gdens_frb.set_xlabel("X [Mpc]", color="white", fontsize=11)
        ax_gdens_frb.set_ylabel("Y [Mpc]", color="white", fontsize=11)
        ax_gdens_frb.tick_params(colors='white', which='both', labelsize=10)  # white numbers and axis hatching
        for spine in ax_gdens_frb.spines.values():
            spine.set_color('#333333')    # dark grey borders
        cbar = fig_gdens_frb.colorbar(im, ax=ax_gdens_frb, fraction=0.046, pad=0.04, shrink=0.8)
        cbar.set_label(r"$\log_{10}$ Gas Density [$g/cm^3$]", color="white", fontsize=10)
        cbar.ax.tick_params(colors='white', labelsize=9)
        cbar.outline.set_color('#333333')
        
        fig_gdens_frb.tight_layout()
        os.makedirs(os.path.join(save_path, "Gas", "GasDensity2D"), exist_ok=True)
        output_name = os.path.join(save_path, "Gas", "GasDensity2D", f"Gas_Density_{i:03d}.png")
        fig_gdens_frb.savefig(output_name, dpi=200, facecolor=fig_gdens_frb.get_facecolor(), bbox_inches="tight")
        plt.close(fig_gdens_frb)
        
        print(f"[{i+1:3d}/{len(plotfiles)}] Saved 2D FRB gas map at z={z_now:.2f}")


        # ---- 2D slice plots ----
        # Global annotations valid for all plots
        slc.annotate_timestamp(corner="upper_left", time=True, draw_inset_box=True)
        slc.annotate_scale(corner="upper_right")
        slc.annotate_grids(alpha=0.3, min_level=1)  # visualize AMR refined grids
        slc._setup_plots() # force matplotlib rendering to manipulate axes and save custom figures

        # Specific settings for each field
        slc.set_cmap(('boxlib', 'density'), cmap="magma")
        # slc.set_zlim(('boxlib', 'density'), 1e-31, 1e-20)  # colorbar CGS range
    
        # Matplotlib logic to insert the particle legend box only on the Density plot
        ax = slc.plots[('boxlib', 'density')].axes
        text_legend = r"$\circ$ DM Particle (CIC)"
        legend_box = AnchoredText(text_legend, loc='lower left', pad=0.7, borderpad=0.7,
                               prop=dict(size=13, color='black'), frameon=True)
        legend_box.patch.set_facecolor('white')
        legend_box.patch.set_alpha(0.95)
        legend_box.patch.set_edgecolor('gray')
        ax.add_artist(legend_box)

    # Create folders (if they do not exist) and save
    os.makedirs(os.path.join(save_path, "Gas", "Density"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "Eint"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "Temp"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "XMom"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "YMom"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "ZMom"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gas", "VelMag"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gravity", "GPot"), exist_ok=True)
    os.makedirs(os.path.join(save_path, "Gravity", "DMmass"), exist_ok=True)

    dens_path = os.path.join(save_path, "Gas", "Density", f"Density_{i:03d}.png")
    eint_path = os.path.join(save_path, "Gas", "Eint", f"InternalEnergy_{i:03d}.png")
    temp_path = os.path.join(save_path, "Gas", "Temp", f"Temperature_{i:03d}.png")
    xmom_path = os.path.join(save_path, "Gas", "XMom", f"Momx_{i:03d}.png")
    ymom_path = os.path.join(save_path, "Gas", "YMom", f"Momy_{i:03d}.png")
    zmom_path = os.path.join(save_path, "Gas", "ZMom", f"Momz_{i:03d}.png")
    velmag_path = os.path.join(save_path, "Gas", "VelMag", f"VelMag_{i:03d}.png")
    gpot_path = os.path.join(save_path, "Gravity", "GPot", f"GPot_{i:03d}.png")
    dmm_path = os.path.join(save_path, "Gravity", "DMmass", f"DMmass_{i:03d}.png")

    slc.plots[('boxlib', 'density')].figure.savefig(dens_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'rho_e')].figure.savefig(eint_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'Temp')].figure.savefig(temp_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'xmom')].figure.savefig(xmom_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'ymom')].figure.savefig(ymom_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'zmom')].figure.savefig(zmom_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'magvel')].figure.savefig(velmag_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'phi_grav')].figure.savefig(gpot_path, dpi=300, bbox_inches='tight')
    slc.plots[('boxlib', 'particle_mass_density')].figure.savefig(dmm_path, dpi=300, bbox_inches='tight')



if __name__ == "__main__":
    main()
    print("All ok.")



# ds.field_list

# [('DM', 'particle_cpu'), ('DM', 'particle_id'), ('DM', 'particle_mass'),
# ('DM', 'particle_position_x'), ('DM', 'particle_position_y'), ('DM', 'particle_position_z'), 
# ('DM', 'particle_xvel'), ('DM', 'particle_yvel'), ('DM', 'particle_zvel'),
#  ('all', 'particle_cpu'), ('all', 'particle_id'), ('all', 'particle_mass'),
#  ('all', 'particle_position_x'), ('all', 'particle_position_y'), 
# ('all', 'particle_position_z'), ('all', 'particle_xvel'), ('all', 'particle_yvel'),
#  ('all', 'particle_zvel'), ('boxlib', 'Ne'), ('boxlib', 'Temp'), ('boxlib', 'density'),
#  ('boxlib', 'grav_x'), ('boxlib', 'grav_y'), ('boxlib', 'grav_z'), ('boxlib', 'magvel'), 
# ('boxlib', 'particle_count'), ('boxlib', 'particle_mass_density'), ('boxlib', 'phi_grav'), 
# ('boxlib', 'pressure'), ('boxlib', 'rho_E'), ('boxlib', 'rho_e'), ('boxlib', 'xmom'), ('boxlib', 'ymom'),
#  ('boxlib', 'zmom'), ('nbody', 'particle_cpu'), ('nbody', 'particle_id'), ('nbody', 'particle_mass'), 
# ('nbody', 'particle_position_x'), ('nbody', 'particle_position_y'), ('nbody', 'particle_position_z'), 
# ('nbody', 'particle_xvel'), ('nbody', 'particle_yvel'), ('nbody', 'particle_zvel')]
