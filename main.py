#!/usr/bin/env python3
"""
Protein Structure Visualization Tool
=====================================
CLI entry point for fetching, analyzing, and visualizing the
SARS-CoV-2 Spike Glycoprotein (PDB: 6VYB) from RCSB PDB.

Usage:
    python main.py                    # Full pipeline with default protein (6VYB)
    python main.py --pdb 6VXX         # Use a different PDB ID
    python main.py --no-3d            # Skip interactive 3D HTML generation
    python main.py --chain A          # Analyze specific chain only
    python main.py --info-only        # Print metadata without generating plots

Author: Naitik Gupta
"""

import argparse
import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.fetcher import fetch_pdb_file, fetch_metadata, fetch_polymer_entities, parse_structure, print_metadata
from src.analyzer import (
    get_chain_stats, get_residue_bfactors, get_residue_composition,
    get_property_distribution, compute_ca_distance_matrix,
    identify_key_regions, summarize_structure
)
from src.visualizer import (
    plot_bfactor_profile, plot_composition, plot_property_distribution,
    plot_distance_matrix, plot_bfactor_by_chain, plot_key_regions_map
)
from src.viewer_3d import create_all_views


def run_pipeline(pdb_id: str = "6VYB", chain_id: str = "A",
                  generate_3d: bool = True, info_only: bool = False):
    """
    Execute the full analysis and visualization pipeline.

    Steps:
        1. Fetch PDB structure and metadata from RCSB
        2. Parse structure with BioPython
        3. Run structural analyses
        4. Generate static visualizations (matplotlib)
        5. Generate interactive 3D viewers (py3Dmol)
    """
    print("=" * 65)
    print("  PROTEIN STRUCTURE VISUALIZATION TOOL")
    print(f"  Target: {pdb_id} | Chain: {chain_id}")
    print("=" * 65)

    # ── Step 1: Fetch ──────────────────────────────────────────
    print("\n[1/5] Fetching structure and metadata...")
    pdb_file = fetch_pdb_file(pdb_id)
    metadata = fetch_metadata(pdb_id)
    print_metadata(metadata)

    try:
        entities = fetch_polymer_entities(pdb_id)
        for ent in entities:
            print(f"  Entity {ent['entity_id']}: {ent['description']}")
            print(f"    Organism: {ent['organism']}")
            print(f"    Sequence length: {ent['sequence_length']} residues")
    except Exception as e:
        print(f"  [warn] Could not fetch entity details: {e}")

    if info_only:
        print("\n[done] Info-only mode — skipping analysis and visualization.")
        return

    # ── Step 2: Parse ──────────────────────────────────────────
    print("\n[2/5] Parsing structure with BioPython...")
    structure = parse_structure(pdb_file, pdb_id)
    model = structure[0]
    print(f"  Chains found: {[c.id for c in model.get_chains()]}")

    # ── Step 3: Analyze ────────────────────────────────────────
    print("\n[3/5] Running structural analyses...")

    chain_stats = get_chain_stats(structure)
    print(f"\n  Chain Statistics:")
    for cs in chain_stats:
        print(f"    Chain {cs['chain_id']}: {cs['residue_count']} residues, "
              f"{cs['atom_count']} atoms, avg B-factor={cs['avg_bfactor']:.1f}")

    bfactors = get_residue_bfactors(structure)
    print(f"  B-factors extracted for {len(bfactors['bfactors'])} residues")

    composition = get_residue_composition(structure)
    print(f"  Amino acid composition: {composition['total_residues']} residues analyzed")

    properties = get_property_distribution(structure)
    for cat, vals in properties["distribution"].items():
        print(f"    {cat}: {vals['count']} ({vals['percentage']}%)")

    key_regions = identify_key_regions(structure, pdb_id)
    if key_regions:
        print(f"\n  Functional Regions ({len(key_regions)} annotated):")
        for r in key_regions:
            print(f"    {r['name']} [{r['start']}-{r['end']}]: {r['description']}")

    # Distance matrix for specified chain
    try:
        dist_data = compute_ca_distance_matrix(structure, chain_id)
        print(f"  Distance matrix: {dist_data['n_residues']}×{dist_data['n_residues']} "
              f"(Chain {chain_id})")
    except Exception as e:
        print(f"  [warn] Could not compute distance matrix for chain {chain_id}: {e}")
        dist_data = None

    # Save analysis summary as JSON
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)
    summary = summarize_structure(structure, pdb_id)
    summary_path = output_dir / f"{pdb_id}_analysis.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n  Analysis summary saved → {summary_path}")

    # ── Step 4: Static Visualizations ──────────────────────────
    print("\n[4/5] Generating static visualizations...")

    plot_bfactor_profile(bfactors, pdb_id, key_regions)
    plot_composition(composition, pdb_id)
    plot_property_distribution(properties, pdb_id)
    plot_bfactor_by_chain(chain_stats, pdb_id)

    if key_regions:
        plot_key_regions_map(key_regions, pdb_id=pdb_id)

    if dist_data:
        plot_distance_matrix(dist_data)

    # ── Step 5: Interactive 3D Views ───────────────────────────
    if generate_3d:
        print("\n[5/5] Generating interactive 3D HTML views...")
        create_all_views(pdb_id)
    else:
        print("\n[5/5] Skipping 3D views (--no-3d flag)")

    print("\n" + "=" * 65)
    print("  PIPELINE COMPLETE")
    print(f"  All outputs saved to: {output_dir.resolve()}")
    print("=" * 65)


def main():
    parser = argparse.ArgumentParser(
        description="Protein Structure Visualization Tool — "
                    "Fetch, analyze, and visualize PDB structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py                       # Analyze SARS-CoV-2 Spike (6VYB)
  python main.py --pdb 6VXX            # Analyze closed-state spike
  python main.py --pdb 7DWZ --chain B  # Different protein, chain B
  python main.py --info-only           # Just print metadata
  python main.py --no-3d               # Skip 3D HTML generation
        """
    )
    parser.add_argument("--pdb", default="6VYB",
                        help="PDB ID to analyze (default: 6VYB)")
    parser.add_argument("--chain", default="A",
                        help="Chain ID for distance matrix (default: A)")
    parser.add_argument("--no-3d", action="store_true",
                        help="Skip interactive 3D HTML generation")
    parser.add_argument("--info-only", action="store_true",
                        help="Only print metadata, skip analysis/visualization")

    args = parser.parse_args()

    run_pipeline(
        pdb_id=args.pdb.upper(),
        chain_id=args.chain,
        generate_3d=not args.no_3d,
        info_only=args.info_only,
    )


if __name__ == "__main__":
    main()
