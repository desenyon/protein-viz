"""
analyzer.py — Structural analysis of protein structures using BioPython.

Computes per-residue B-factors, secondary structure distribution,
chain statistics, residue composition, and distance matrices for
detailed structural characterization.
"""

import numpy as np
from collections import Counter, defaultdict
from typing import Optional
from Bio.PDB import Structure, PPBuilder, is_aa
from Bio.PDB.DSSP import DSSP


# Standard amino acid 3-letter to 1-letter mapping
AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Amino acid property classifications
HYDROPHOBIC = set("AILMFWVP")
POLAR = set("STNQY")
CHARGED_POS = set("RHK")
CHARGED_NEG = set("DE")
SPECIAL = set("CGP")


def get_chain_stats(structure: Structure) -> list:
    """
    Compute per-chain statistics: residue count, atom count,
    average B-factor, and sequence.
    """
    model = structure[0]
    stats = []

    for chain in model:
        residues = [r for r in chain.get_residues() if is_aa(r, standard=True)]
        atoms = list(chain.get_atoms())
        b_factors = [a.get_bfactor() for a in atoms]

        seq = "".join(
            AA_MAP.get(r.get_resname(), "X") for r in residues
        )

        stats.append({
            "chain_id": chain.id,
            "residue_count": len(residues),
            "atom_count": len(atoms),
            "avg_bfactor": round(np.mean(b_factors), 2) if b_factors else 0,
            "min_bfactor": round(np.min(b_factors), 2) if b_factors else 0,
            "max_bfactor": round(np.max(b_factors), 2) if b_factors else 0,
            "sequence_length": len(seq),
        })

    return stats


def get_residue_bfactors(structure: Structure, chain_id: Optional[str] = None) -> dict:
    """
    Extract per-residue average B-factors.

    Returns dict with:
      - residue_ids: list of (chain_id, resseq) tuples
      - residue_names: list of 3-letter residue names
      - bfactors: list of average B-factors per residue
    """
    model = structure[0]
    result = {"residue_ids": [], "residue_names": [], "bfactors": []}

    chains = [model[chain_id]] if chain_id else model.get_chains()

    for chain in chains:
        for residue in chain.get_residues():
            if not is_aa(residue, standard=True):
                continue

            atoms = list(residue.get_atoms())
            avg_bf = np.mean([a.get_bfactor() for a in atoms])

            result["residue_ids"].append((chain.id, residue.id[1]))
            result["residue_names"].append(residue.get_resname())
            result["bfactors"].append(round(avg_bf, 2))

    return result


def get_residue_composition(structure: Structure) -> dict:
    """
    Compute amino acid composition across all chains.

    Returns dict mapping 1-letter codes to counts and percentages.
    """
    model = structure[0]
    counter = Counter()

    for chain in model:
        for residue in chain.get_residues():
            if is_aa(residue, standard=True):
                one_letter = AA_MAP.get(residue.get_resname(), "X")
                counter[one_letter] += 1

    total = sum(counter.values())
    composition = {}
    for aa, count in sorted(counter.items()):
        composition[aa] = {
            "count": count,
            "percentage": round(100 * count / total, 2) if total else 0,
        }

    return {"composition": composition, "total_residues": total}


def get_property_distribution(structure: Structure) -> dict:
    """
    Classify residues by biochemical property:
    hydrophobic, polar, positively charged, negatively charged, special.
    """
    model = structure[0]
    categories = {
        "Hydrophobic": 0,
        "Polar": 0,
        "Positive (+)": 0,
        "Negative (-)": 0,
        "Special": 0,
    }

    total = 0
    for chain in model:
        for residue in chain.get_residues():
            if not is_aa(residue, standard=True):
                continue
            aa = AA_MAP.get(residue.get_resname(), "X")
            total += 1
            if aa in HYDROPHOBIC:
                categories["Hydrophobic"] += 1
            elif aa in POLAR:
                categories["Polar"] += 1
            elif aa in CHARGED_POS:
                categories["Positive (+)"] += 1
            elif aa in CHARGED_NEG:
                categories["Negative (-)"] += 1
            else:
                categories["Special"] += 1

    # Convert to percentages
    distribution = {}
    for cat, count in categories.items():
        distribution[cat] = {
            "count": count,
            "percentage": round(100 * count / total, 2) if total else 0,
        }

    return {"distribution": distribution, "total_residues": total}


def compute_ca_distance_matrix(structure: Structure, chain_id: str = "A",
                                max_residues: int = 200) -> dict:
    """
    Compute Cα-Cα distance matrix for a chain.

    Returns dict with distance matrix (numpy array) and residue labels.
    Limits to first `max_residues` for performance.
    """
    model = structure[0]
    chain = model[chain_id]

    ca_atoms = []
    labels = []

    for residue in chain.get_residues():
        if not is_aa(residue, standard=True):
            continue
        if "CA" in residue:
            ca_atoms.append(residue["CA"].get_vector().get_array())
            resname = AA_MAP.get(residue.get_resname(), "X")
            labels.append(f"{resname}{residue.id[1]}")

        if len(ca_atoms) >= max_residues:
            break

    coords = np.array(ca_atoms)
    n = len(coords)

    # Vectorized distance computation
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))

    return {
        "distance_matrix": dist_matrix,
        "labels": labels,
        "n_residues": n,
        "chain_id": chain_id,
    }


def identify_key_regions(structure: Structure, pdb_id: str = "6VYB") -> list:
    """
    Return annotated functional regions for well-known proteins.
    Currently supports SARS-CoV-2 Spike (6VYB).
    """
    regions = {
        "6VYB": [
            {"name": "N-terminal Domain (NTD)", "start": 14, "end": 305,
             "description": "Involved in host cell attachment and immune evasion"},
            {"name": "Receptor Binding Domain (RBD)", "start": 319, "end": 541,
             "description": "Directly binds ACE2 receptor; primary target for neutralizing antibodies"},
            {"name": "Receptor Binding Motif (RBM)", "start": 437, "end": 508,
             "description": "Contact interface with ACE2; contains key mutation sites"},
            {"name": "Furin Cleavage Site (S1/S2)", "start": 682, "end": 685,
             "description": "RRAR motif; unique to SARS-CoV-2, enhances cell entry"},
            {"name": "Fusion Peptide (FP)", "start": 788, "end": 806,
             "description": "Inserts into host cell membrane during fusion"},
            {"name": "Heptad Repeat 1 (HR1)", "start": 912, "end": 984,
             "description": "Forms coiled-coil structure critical for membrane fusion"},
            {"name": "Central Helix (CH)", "start": 986, "end": 1035,
             "description": "Connects HR1 to connector domain"},
        ],
    }

    return regions.get(pdb_id.upper(), [])


def summarize_structure(structure: Structure, pdb_id: str = "6VYB") -> dict:
    """
    Generate a full structural summary combining all analyses.
    """
    chain_stats = get_chain_stats(structure)
    composition = get_residue_composition(structure)
    properties = get_property_distribution(structure)
    regions = identify_key_regions(structure, pdb_id)

    return {
        "pdb_id": pdb_id,
        "chain_stats": chain_stats,
        "composition": composition,
        "property_distribution": properties,
        "key_regions": regions,
        "total_chains": len(chain_stats),
        "total_atoms": sum(c["atom_count"] for c in chain_stats),
        "total_residues": composition["total_residues"],
    }
