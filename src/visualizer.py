"""
visualizer.py — Static 2D/3D visualizations of protein structures.

Generates publication-quality matplotlib figures:
  - B-factor heatmaps
  - Amino acid composition bar charts
  - Property distribution pie charts
  - Cα distance matrices (contact maps)
  - Ramachandran-like residue position plots

All outputs are saved to the output/ directory.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from typing import Optional

OUTPUT_DIR = Path(__file__).parent.parent / "output"


def setup_style():
    """Apply a clean, publication-quality plot style."""
    plt.rcParams.update({
        "figure.facecolor": "#0d1117",
        "axes.facecolor": "#161b22",
        "axes.edgecolor": "#30363d",
        "axes.labelcolor": "#c9d1d9",
        "text.color": "#c9d1d9",
        "xtick.color": "#8b949e",
        "ytick.color": "#8b949e",
        "grid.color": "#21262d",
        "font.family": "monospace",
        "font.size": 11,
    })


def plot_bfactor_profile(bfactor_data: dict, pdb_id: str = "6VYB",
                          key_regions: list = None,
                          filename: str = "bfactor_profile.png") -> Path:
    """
    Plot per-residue B-factor profile with optional functional region annotations.

    B-factors indicate atomic displacement / flexibility.
    High B-factors → flexible/disordered regions.
    Low B-factors → rigid/well-ordered regions.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    bfactors = np.array(bfactor_data["bfactors"])
    residue_ids = bfactor_data["residue_ids"]
    x = np.arange(len(bfactors))

    fig, ax = plt.subplots(figsize=(16, 6))

    # Color by B-factor magnitude
    colors = plt.cm.coolwarm(bfactors / max(bfactors.max(), 1))
    ax.bar(x, bfactors, color=colors, width=1.0, edgecolor="none", alpha=0.85)

    # Annotate key regions if provided
    region_colors = ["#58a6ff", "#f0883e", "#d2a8ff", "#3fb950",
                     "#f778ba", "#79c0ff", "#ffa657"]
    if key_regions:
        for i, region in enumerate(key_regions):
            color = region_colors[i % len(region_colors)]
            # Find residue indices in our data that fall in this region
            start_idx = None
            end_idx = None
            for j, (chain_id, resseq) in enumerate(residue_ids):
                if resseq >= region["start"] and start_idx is None:
                    start_idx = j
                if resseq <= region["end"]:
                    end_idx = j

            if start_idx is not None and end_idx is not None:
                ax.axvspan(start_idx, end_idx, alpha=0.15, color=color,
                          label=region["name"])

    ax.set_xlabel("Residue Index", fontsize=12)
    ax.set_ylabel("B-factor (Å²)", fontsize=12)
    ax.set_title(f"Per-Residue B-factor Profile — {pdb_id} (SARS-CoV-2 Spike)",
                 fontsize=14, fontweight="bold", pad=15)

    if key_regions:
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95,
                 facecolor="#161b22", edgecolor="#30363d",
                 bbox_to_anchor=(0.0, 1.0))

    ax.set_xlim(0, len(bfactors))
    fig.subplots_adjust(top=0.88)
    fig.tight_layout()

    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved B-factor profile → {path}")
    return path


def plot_composition(composition_data: dict, pdb_id: str = "6VYB",
                     filename: str = "aa_composition.png") -> Path:
    """
    Bar chart of amino acid composition with color-coded properties.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    comp = composition_data["composition"]
    aas = sorted(comp.keys())
    counts = [comp[aa]["count"] for aa in aas]
    pcts = [comp[aa]["percentage"] for aa in aas]

    # Color by property
    hydrophobic = set("AILMFWVP")
    polar = set("STNQY")
    positive = set("RHK")
    negative = set("DE")

    colors = []
    for aa in aas:
        if aa in hydrophobic:
            colors.append("#f0883e")   # orange
        elif aa in polar:
            colors.append("#58a6ff")   # blue
        elif aa in positive:
            colors.append("#3fb950")   # green
        elif aa in negative:
            colors.append("#f85149")   # red
        else:
            colors.append("#8b949e")   # gray

    fig, ax = plt.subplots(figsize=(14, 6))
    bars = ax.bar(aas, counts, color=colors, edgecolor="#30363d", linewidth=0.5)

    # Add percentage labels on top
    for bar, pct in zip(bars, pcts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                f"{pct:.1f}%", ha="center", va="bottom", fontsize=8,
                color="#8b949e")

    # Legend
    legend_patches = [
        mpatches.Patch(color="#f0883e", label="Hydrophobic"),
        mpatches.Patch(color="#58a6ff", label="Polar"),
        mpatches.Patch(color="#3fb950", label="Positive (+)"),
        mpatches.Patch(color="#f85149", label="Negative (-)"),
        mpatches.Patch(color="#8b949e", label="Special (G, C)"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9,
             framealpha=0.8, facecolor="#161b22", edgecolor="#30363d")

    ax.set_xlabel("Amino Acid", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(f"Amino Acid Composition — {pdb_id} (SARS-CoV-2 Spike)",
                 fontsize=14, fontweight="bold", pad=15)

    fig.tight_layout()
    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved composition chart → {path}")
    return path


def plot_property_distribution(property_data: dict, pdb_id: str = "6VYB",
                                filename: str = "property_dist.png") -> Path:
    """
    Donut chart showing biochemical property distribution of residues.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    dist = property_data["distribution"]
    labels = list(dist.keys())
    sizes = [dist[l]["count"] for l in labels]
    pcts = [dist[l]["percentage"] for l in labels]

    colors = ["#f0883e", "#58a6ff", "#3fb950", "#f85149", "#8b949e"]

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=None, autopct=lambda p: f"{p:.1f}%",
        colors=colors, startangle=90, pctdistance=0.8,
        wedgeprops=dict(width=0.4, edgecolor="#0d1117", linewidth=2)
    )

    for t in autotexts:
        t.set_fontsize(10)
        t.set_color("#c9d1d9")

    # Center label
    ax.text(0, 0, f"{pdb_id}\n{property_data['total_residues']}\nresidues",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color="#c9d1d9")

    ax.legend(labels, loc="lower right", fontsize=10,
             framealpha=0.8, facecolor="#161b22", edgecolor="#30363d")
    ax.set_title(f"Residue Property Distribution — {pdb_id}",
                 fontsize=14, fontweight="bold", pad=15)

    fig.tight_layout()
    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved property distribution → {path}")
    return path


def plot_distance_matrix(dist_data: dict, filename: str = "contact_map.png") -> Path:
    """
    Heatmap of Cα-Cα distances (contact map).

    Residue pairs < 8Å apart are considered in contact.
    Reveals domain boundaries and structural organization.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    matrix = dist_data["distance_matrix"]
    chain_id = dist_data["chain_id"]
    n = dist_data["n_residues"]

    fig, ax = plt.subplots(figsize=(10, 10))

    # Use a custom colormap: close contacts in bright colors
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=50, aspect="equal",
                   interpolation="nearest")

    # Add contact threshold line (8 Å)
    contact_map = (matrix < 8.0).astype(float)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="Distance (Å)")
    cbar.ax.yaxis.label.set_color("#c9d1d9")
    cbar.ax.tick_params(colors="#8b949e")

    # Label every 50th residue
    tick_step = max(1, n // 10)
    tick_positions = range(0, n, tick_step)
    tick_labels = [dist_data["labels"][i] for i in tick_positions]
    ax.set_xticks(list(tick_positions))
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(list(tick_positions))
    ax.set_yticklabels(tick_labels, fontsize=7)

    ax.set_xlabel("Residue", fontsize=12)
    ax.set_ylabel("Residue", fontsize=12)
    ax.set_title(f"Cα Distance Matrix — Chain {chain_id} (first {n} residues)",
                 fontsize=14, fontweight="bold", pad=15)

    fig.tight_layout()
    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved distance matrix → {path}")
    return path


def plot_bfactor_by_chain(chain_stats: list, pdb_id: str = "6VYB",
                           filename: str = "chain_bfactors.png") -> Path:
    """
    Comparative B-factor statistics across chains.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    chains = [c["chain_id"] for c in chain_stats]
    avgs = [c["avg_bfactor"] for c in chain_stats]
    mins = [c["min_bfactor"] for c in chain_stats]
    maxs = [c["max_bfactor"] for c in chain_stats]

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(chains))
    width = 0.25

    ax.bar(x - width, mins, width, label="Min", color="#3fb950", alpha=0.8)
    ax.bar(x, avgs, width, label="Average", color="#58a6ff", alpha=0.8)
    ax.bar(x + width, maxs, width, label="Max", color="#f85149", alpha=0.8)

    ax.set_xlabel("Chain", fontsize=12)
    ax.set_ylabel("B-factor (Å²)", fontsize=12)
    ax.set_title(f"B-factor Statistics by Chain — {pdb_id}",
                 fontsize=14, fontweight="bold", pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(chains, fontsize=12)
    ax.legend(framealpha=0.8, facecolor="#161b22", edgecolor="#30363d")

    fig.tight_layout()
    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved chain B-factors → {path}")
    return path


def plot_key_regions_map(key_regions: list, total_length: int = 1273,
                          pdb_id: str = "6VYB",
                          filename: str = "domain_map.png") -> Path:
    """
    Linear domain map showing functional regions along the sequence.
    """
    setup_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16, 5))

    # Draw full-length bar
    ax.barh(0, total_length, height=0.5, color="#21262d", edgecolor="#30363d")

    region_colors = ["#79c0ff", "#ffa657", "#d2a8ff", "#7ee787",
                     "#ff9bce", "#a5d6ff", "#ffc680"]

    for i, region in enumerate(key_regions):
        color = region_colors[i % len(region_colors)]
        start = region["start"]
        width = region["end"] - region["start"]

        ax.barh(0, width, left=start, height=0.5, color=color,
               edgecolor="#0d1117", linewidth=1, alpha=0.9)

    # Place labels below the bar with vertical offset to avoid overlap
    label_y_offsets = [-0.55, -0.75, -0.95, -0.55, -0.75, -0.55, -0.75]
    for i, region in enumerate(key_regions):
        color = region_colors[i % len(region_colors)]
        start = region["start"]
        width = region["end"] - region["start"]
        mid = start + width / 2
        y_off = label_y_offsets[i % len(label_y_offsets)]

        # Draw connector line
        ax.plot([mid, mid], [-0.25, y_off + 0.08], color=color,
               linewidth=0.8, alpha=0.6)
        ax.text(mid, y_off, region["name"], ha="center", va="top",
               fontsize=8, color=color, fontweight="bold")

    ax.set_xlim(0, total_length)
    ax.set_ylim(-1.3, 0.8)
    ax.set_xlabel("Residue Position", fontsize=12)
    ax.set_title(f"Functional Domain Map — {pdb_id} (SARS-CoV-2 Spike Glycoprotein)",
                 fontsize=14, fontweight="bold", pad=15)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)

    fig.tight_layout()
    path = OUTPUT_DIR / filename
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] Saved domain map → {path}")
    return path
