# audiogram_figures.py
# Author: Narges Shahmohammad
# Affiliation: Cellular & Molecular Research Center,
#   Research Institute for Prevention of Non-Communicable Diseases,
#   Qazvin University of Medical Sciences, Qazvin, Iran
# Contact: n.shahmohammad@qums.ac.ir
#
# Purpose:
#   Generate publication-ready pure-tone audiogram figures:
#     1) Single-subject audiograms (Variant B: lines + large O/X markers)
#     2) A two-panel audiogram for siblings (II-2 & II-3)
#     3) A 1×2 composite: (a) Pedigree + (b) II-1 (Proband) audiogram
#
# Key visual choices (per manuscript style):
#   - Severity bands (Normal…Profound) are shaded (light gray) with labels at the right edge.
#   - Legend (O/X) is placed at the upper-left.
#   - Titles:
#       • Singles: use the pedigree number only (e.g., "II-1", "II-2", "II-3")
#       • Siblings two-panel: each panel title is just the pedigree number; panel letters (a), (b) optional
#       • Pedigree+Proband composite: panel titles "(a) Pedigree" and "(b) II-1";
#         the internal title inside the Proband audiogram is removed to avoid duplication.
#
# Python: 3.11+
# Dependencies: matplotlib, numpy
#
# License:
#   This script is released for academic and research purposes only.
#   Please cite the related article if you use or modify this code.

from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg  # only for composing with pedigree

# ----------------------------
# Config (edit as needed)
# ----------------------------

# Octave frequencies (Hz) shown at equal spacing (not logarithmic)
FREQ_LABELS = ['250', '500', '1000', '2000', '4000', '8000']
X = np.arange(len(FREQ_LABELS))  # fixed positions for octave spacing

# Audiometric thresholds (example data; replace with your measurements)
# Order per ear: 250, 500, 1000, 2000, 4000, 8000 Hz  (dB HL)
SUBJECTS = {
    "II-1": {  # Proband
        "right": [40, 45, 50, 55, 55, 40],
        "left":  [40, 50, 50, 55, 60, 45],
    },
    "II-2": {
        "right": [40, 50, 55, 60, 60, 65],
        "left":  [45, 50, 55, 55, 60, 65],
    },
    "II-3": {
        "right": [25, 30, 45, 50, 50, 50],
        "left":  [35, 40, 55, 55, 55, 55],
    },
}

# Severity bands (dB HL) for shading + right-side labels
# Display text can be localized/changed if needed.
SEVERITY_BANDS = [
    ("Normal",   -10, 20),
    ("Mild",       21, 40),
    ("Moderate",   41, 70),
    ("Severe",     71, 95),
    ("Profound",   96, 120),  # Profound (≥96 dB HL)
]

# Paths
OUTDIR = Path("fig_out_final")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Provide your pedigree image path (PNG/TIFF). If not present, the 1×2 composite is skipped.
PEDIGREE_IMAGE_PATH = Path("pedigree.png")

# Figure output formats
EXPORT_EXTS = ("png", "tiff", "pdf", "svg")

# Styling (global)
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.facecolor": "white",
    "grid.color": "0.75",
    "grid.linestyle": "--",
    "grid.linewidth": 0.4,
})

# ----------------------------
# Helpers
# ----------------------------

def add_severity_bands_with_right_labels(ax: plt.Axes) -> None:
    """
    Shade severity bands and place their labels at the right side, centered in each band.
    """
    for i, (name, y0, y1) in enumerate(SEVERITY_BANDS):
        # alternating light gray for print visibility
        ax.axhspan(y0, y1, color=('0.92' if i % 2 == 0 else '0.97'), zorder=0)
        ax.text(len(X) + 0.12, (y0 + y1) / 2, name,
                va='center', fontsize=8, color='0.25')

def save_all(fig: plt.Figure, stub: str | Path) -> None:
    """
    Save figure to OUTDIR in multiple formats.
    """
    stub = OUTDIR / str(stub)
    for ext in EXPORT_EXTS:
        fig.savefig(f"{stub}.{ext}", bbox_inches='tight')

def plot_single_variantB(number: str, right: list[float], left: list[float],
                         show_title_inside: bool = True,
                         stub: str | Path = "audiogram_single") -> None:
    """
    Single audiogram (Variant B: lines + large markers, O/X).
    - Shading + right-side labels
    - Legend at upper-left
    - Title: pedigree number only (if show_title_inside=True)
    """
    fig, ax = plt.subplots(figsize=(5.8, 4.8), dpi=300)
    add_severity_bands_with_right_labels(ax)

    # Grid, axes
    ax.set_axisbelow(True)
    ax.grid(True, linestyle='--', linewidth=0.4)

    # Curves (Right=O, Left=X) — bigger markers so O/X pop above shading
    ax.plot(X, right, 'o-', label='Right (O)', linewidth=1.8, markersize=9, markeredgewidth=1.6, zorder=3)
    ax.plot(X, left,  'x--', label='Left (X)',  linewidth=1.5, markersize=10, markeredgewidth=1.8, zorder=3)

    # Axes labeling
    ax.set_xticks(X); ax.set_xticklabels(FREQ_LABELS)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Hearing Level (dB HL)")
    ax.set_ylim(120, -10)   # invert dB HL
    ax.set_xlim(-0.2, len(X) - 0.2)

    # Legend and (optional) title
    ax.legend(loc='upper left', frameon=True)
    if show_title_inside:
        ax.set_title(number, weight='bold')

    fig.tight_layout()
    save_all(fig, stub)
    plt.close(fig)

def plot_two_panel_variantB(names: list[str], letters: list[str] | None,
                            stub: str | Path) -> None:
    """
    Two-panel audiogram (e.g., siblings II-2 & II-3).
    - Each panel titled with the pedigree number (e.g., "II-2")
    - Optionally prefix titles with panel letters (a), (b)
    - Legends on both panels (upper-left)
    """
    assert len(names) == 2, "Provide exactly two pedigree numbers"
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), dpi=300, sharey=True)

    for idx, (ax, num) in enumerate(zip(axes, names)):
        d = SUBJECTS[num]
        add_severity_bands_with_right_labels(ax)
        ax.set_axisbelow(True)
        ax.grid(True, linestyle='--', linewidth=0.4)
        ax.plot(X, d["right"], 'o-', label='Right (O)', linewidth=1.8, markersize=9, markeredgewidth=1.6, zorder=3)
        ax.plot(X, d["left"],  'x--', label='Left (X)',  linewidth=1.5, markersize=10, markeredgewidth=1.8, zorder=3)
        ax.set_xticks(X); ax.set_xticklabels(FREQ_LABELS)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylim(120, -10); ax.set_xlim(-0.2, len(X)-0.2)
        # Title
        title = num if not letters else f"({letters[idx]}) {num}"
        ax.set_title(title, weight='bold')
        # Legend
        ax.legend(loc='upper left', frameon=True)

    axes[0].set_ylabel("Hearing Level (dB HL)")
    fig.tight_layout()
    save_all(fig, stub)
    plt.close(fig)

def compose_pedigree_with_proband(pedigree_path: Path,
                                  proband_stub_png: Path,
                                  out_stub: str | Path,
                                  letters: tuple[str, str] = ("a", "b")) -> None:
    """
    Compose a 1×2 figure: (a) Pedigree image + (b) Proband audiogram image (PNG).
    Assumes the Proband audiogram was saved WITHOUT an internal title to avoid duplication.
    Panel titles are set as "(a) Pedigree" and "(b) II-1".
    """
    if not (pedigree_path.exists() and proband_stub_png.exists()):
        print("[INFO] Skipping composite: missing image(s).")
        return

    ped_img = mpimg.imread(str(pedigree_path))
    aud_img = mpimg.imread(str(proband_stub_png))

    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=300)
    axs[0].imshow(ped_img); axs[0].axis('off'); axs[0].set_title(f"({letters[0]}) Pedigree", weight='bold')
    axs[1].imshow(aud_img); axs[1].axis('off'); axs[1].set_title(f"({letters[1]}) II-1", weight='bold')
    fig.tight_layout()
    save_all(fig, out_stub)
    plt.close(fig)

# ----------------------------
# Main (runs all outputs)
# ----------------------------
if __name__ == "__main__":
    # 1) Singles — title = pedigree number only
    for num, d in SUBJECTS.items():
        stub = f"audiogram_{num}_variantB_titleNumber"
        plot_single_variantB(num, d["right"], d["left"], show_title_inside=True, stub=stub)

    # 2) Two-panel siblings (II-2 & II-3) — titles = numbers, add panel letters (a), (b)
    plot_two_panel_variantB(names=["II-2", "II-3"],
                            letters=["a", "b"],
                            stub="audiograms_siblings_II2_II3_variantB_titleNumbers_withLetters")

    # 3) Proband audiogram without internal title (to avoid duplicate title in composite)
    plot_single_variantB("II-1",
                         SUBJECTS["II-1"]["right"], SUBJECTS["II-1"]["left"],
                         show_title_inside=False,
                         stub="audiogram_II-1_variantB_noInnerTitle")

    # 4) 1×2 composite: (a) Pedigree + (b) II-1 audiogram
    compose_pedigree_with_proband(
        pedigree_path=PEDIGREE_IMAGE_PATH,
        proband_stub_png=OUTDIR / "audiogram_II-1_variantB_noInnerTitle.png",
        out_stub="Figure_proband_with_pedigree_1x2_withLetters",
        letters=("a", "b"),
    )

    print(f"Done. Figures saved to: {OUTDIR.resolve()}")
