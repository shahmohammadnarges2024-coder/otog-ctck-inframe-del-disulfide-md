# audiogram_plotter.py
# Author: Narges Shahmohammad
# Affiliation: Cellular and Molecular Research Center,
#   Research Institute for Prevention of Non-Communicable Diseases,
#   Qazvin University of Medical Sciences, Qazvin, Iran
# Contact: n.shahmohammad@qums.ac.ir
#
# Purpose: Redraw pure-tone audiograms for publication-quality figures
# Version: 1.0
# Date: September 2025
# Python: 3.11.8
# Dependencies: Matplotlib 3.7.2, NumPy
#
# License: This script is released for academic and research purposes only.
#          Please cite the related article if you use or modify this code.

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg  # only needed for pedigree combination

# ----------------------------
# Configurable inputs
# ----------------------------

# Octave frequencies (Hz) shown at equal spacing (not logarithmic)
FREQ_LABELS = ['250','500','1000','2000','4000','8000']
X = np.arange(len(FREQ_LABELS))  # fixed positions for octave spacing

# Example data (replace with your measurements)
SUBJECTS = {
    "Proband (II-1)": {
        "right": [40, 45, 50, 55, 55, 40],
        "left":  [40, 50, 50, 55, 60, 45]
    },
    "II-2": {
        "right": [40, 50, 55, 60, 60, 65],
        "left":  [45, 50, 55, 55, 60, 65]
    },
    "II-3": {
        "right": [25, 30, 45, 50, 50, 50],
        "left":  [35, 40, 55, 55, 55, 55]
    }
}

# Severity bands to shade (dB HL)
SEVERITY_BANDS = [
    ("Normal",   -10, 20),
    ("Mild",       21, 40),
    ("Moderate",   41, 70),
    ("Severe",     71, 95),
    ("Profound",   96, 120),
]

# Output directory
OUTDIR = Path("fig_out")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Optional pedigree image (PNG/TIFF) for 2x2 combined figure
PEDIGREE_IMAGE_PATH = "pedigree.png"   # put your pedigree image here (PNG/TIFF)

# ----------------------------
# Plotting helpers
# ----------------------------

def add_severity_bands(ax):
    """Shade horizontal severity bands and label at right side."""
    for name, y0, y1 in SEVERITY_BANDS:
        ax.axhspan(y0, y1, color='0.92', zorder=0)
        # Put labels slightly outside the right edge
        ax.text(len(X) + 0.12, (y0 + y1) / 2, name, va='center', fontsize=8)

def plot_single_audiogram(title, right, left, save_stub):
    """Draw one subject's audiogram and save as PNG/TIFF (600 dpi)."""
    fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=300)
    add_severity_bands(ax)

    # Grid (publication-friendly)
    ax.set_axisbelow(True)
    ax.grid(True, which='both', linestyle='--', linewidth=0.4)

    # Threshold curves (O = right, X = left)
    ax.plot(X, right, 'o-', label='Right (O)', linewidth=1.5, markersize=5)
    ax.plot(X, left,  'x--', label='Left (X)',  linewidth=1.2, markersize=6)

    # Axes and labels
    ax.set_xticks(X)
    ax.set_xticklabels(FREQ_LABELS)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Hearing Level (dB HL)")
    ax.set_ylim(120, -10)                     # invert dB HL axis
    ax.set_xlim(-0.2, len(X) - 0.2)
    ax.legend(loc='lower right', fontsize=8)
    ax.set_title(f"Pure-tone Audiogram – {title}", fontsize=10, weight='bold')

    fig.tight_layout()

    png = OUTDIR / f"{save_stub}.png"
    tiff = OUTDIR / f"{save_stub}.tiff"
    fig.savefig(png,  dpi=600, bbox_inches='tight')
    fig.savefig(tiff, dpi=600, bbox_inches='tight')
    plt.close(fig)
    return png, tiff

def plot_three_panel(subject_dict, save_stub):
    """Three subjects in one row (for supplementary)."""
    names = list(subject_dict.keys())
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), dpi=300, sharey=True)

    for ax, name in zip(axes, names):
        add_severity_bands(ax)
        ax.grid(True, linestyle='--', linewidth=0.4)
        ax.plot(X, subject_dict[name]["right"], 'o-', linewidth=1.5, markersize=5, label='Right (O)')
        ax.plot(X, subject_dict[name]["left"],  'x--', linewidth=1.2, markersize=6, label='Left (X)')
        ax.set_xticks(X)
        ax.set_xticklabels(FREQ_LABELS)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylim(120, -10)
        ax.set_xlim(-0.2, len(X) - 0.2)
        ax.set_title(name, fontsize=10, weight='bold')

    axes[0].set_ylabel("Hearing Level (dB HL)")
    axes[-1].legend(loc='lower right', fontsize=8)
    fig.tight_layout()

    png = OUTDIR / f"{save_stub}.png"
    tiff = OUTDIR / f"{save_stub}.tiff"
    fig.savefig(png,  dpi=600, bbox_inches='tight')
    fig.savefig(tiff, dpi=600, bbox_inches='tight')
    plt.close(fig)
    return png, tiff

def plot_2x2_with_pedigree(pedigree_path, subject_dict, save_stub):
    """2×2 composite: (a) pedigree + (b,c,d) three audiograms."""
    # Ensure individual panels exist (or create inline)
    panels = []
    for name in subject_dict:
        stub = name.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")
        png, _ = plot_single_audiogram(name, subject_dict[name]["right"], subject_dict[name]["left"], f"audiogram_{stub}")
        panels.append(png)

    # Load images
    ped_img = mpimg.imread(pedigree_path)
    imgs = [mpimg.imread(str(p)) for p in panels]

    fig, axs = plt.subplots(2, 2, figsize=(12, 10), dpi=300)
    axs[0,0].imshow(ped_img); axs[0,0].axis('off'); axs[0,0].set_title("(a) Pedigree", fontsize=11, weight='bold')
    axs[0,1].imshow(imgs[0]);  axs[0,1].axis('off'); axs[0,1].set_title("(b) " + list(subject_dict.keys())[0], fontsize=11, weight='bold')
    axs[1,0].imshow(imgs[1]);  axs[1,0].axis('off'); axs[1,0].set_title("(c) " + list(subject_dict.keys())[1], fontsize=11, weight='bold')
    axs[1,1].imshow(imgs[2]);  axs[1,1].axis('off'); axs[1,1].set_title("(d) " + list(subject_dict.keys())[2], fontsize=11, weight='bold')

    fig.tight_layout()
    png = OUTDIR / f"{save_stub}.png"
    tiff = OUTDIR / f"{save_stub}.tiff"
    fig.savefig(png,  dpi=600, bbox_inches='tight')
    fig.savefig(tiff, dpi=600, bbox_inches='tight')
    plt.close(fig)
    return png, tiff

# ----------------------------
# Main run (example)
# ----------------------------
if __name__ == "__main__":
    # 1) Save individual audiograms
    for name, data in SUBJECTS.items():
        stub = name.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")
        plot_single_audiogram(name, data["right"], data["left"], f"audiogram_{stub}")

    # 2) Save a three-panel combined audiogram (optional)
    plot_three_panel(SUBJECTS, "audiograms_three_panel")

    # 3) Save a 2×2 composite with pedigree (requires PEDIGREE_IMAGE_PATH to exist)
    if os.path.exists(PEDIGREE_IMAGE_PATH):
        plot_2x2_with_pedigree(PEDIGREE_IMAGE_PATH, SUBJECTS, "Figure1_pedigree_audiograms_2x2")
    else:
        print(f"[INFO] Skipping 2×2 composite: pedigree image not found at {PEDIGREE_IMAGE_PATH}")
