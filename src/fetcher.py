"""
fetcher.py — Fetch and cache PDB structures from RCSB Protein Data Bank.

Supports downloading PDB files and querying the RCSB REST API for
metadata (resolution, method, organism, citation, etc.).
"""

import os
import json
import requests
from pathlib import Path
from typing import Optional
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


DATA_DIR = Path(__file__).parent.parent / "data"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"
RCSB_API_URL = "https://data.rcsb.org/rest/v1/core"


def fetch_pdb_file(pdb_id: str, file_format: str = "pdb", force: bool = False) -> Path:
    """
    Download a structure file from RCSB PDB.

    Args:
        pdb_id: 4-character PDB identifier (e.g., '6VYB').
        file_format: 'pdb' or 'cif' (mmCIF).
        force: If True, re-download even if cached.

    Returns:
        Path to the downloaded file.
    """
    pdb_id = pdb_id.upper()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    ext = "pdb" if file_format == "pdb" else "cif"
    local_path = DATA_DIR / f"{pdb_id}.{ext}"

    if local_path.exists() and not force:
        print(f"[cache] Using cached {local_path.name}")
        return local_path

    url = f"{RCSB_DOWNLOAD_URL}/{pdb_id}.{ext}"
    print(f"[fetch] Downloading {pdb_id}.{ext} from RCSB PDB...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    local_path.write_text(resp.text)
    print(f"[fetch] Saved to {local_path}")
    return local_path


def fetch_metadata(pdb_id: str) -> dict:
    """
    Query the RCSB REST API for entry-level metadata.

    Returns a dict with title, resolution, method, deposition date,
    authors, organism, molecular weight, keywords, and citation info.
    """
    pdb_id = pdb_id.upper()
    url = f"{RCSB_API_URL}/entry/{pdb_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    # Extract structured metadata
    struct = data.get("struct", {})
    info = data.get("rcsb_entry_info", {})
    citation = data.get("rcsb_primary_citation", {})
    em = data.get("em3d_reconstruction", [{}])
    em_info = em[0] if isinstance(em, list) and em else em

    meta = {
        "pdb_id": pdb_id,
        "title": struct.get("title", "N/A"),
        "resolution": info.get("resolution_combined", ["N/A"]),
        "experimental_method": info.get("experimental_method", "N/A"),
        "deposition_date": data.get("rcsb_accession_info", {}).get("deposit_date", "N/A"),
        "molecular_weight_kda": round(info.get("molecular_weight", 0), 2),
        "polymer_entity_count": info.get("polymer_entity_count", 0),
        "deposited_atom_count": info.get("deposited_atom_count", 0),
        "deposited_residue_count": info.get("deposited_polymer_monomer_count", 0),
        "disulfide_bonds": info.get("disulfide_bond_count", 0),
        "keywords": data.get("struct_keywords", {}).get("text", "N/A"),
        "authors": [a.get("name", "") for a in data.get("audit_author", [])],
        "citation_title": citation.get("title", "N/A"),
        "citation_journal": citation.get("journal_abbrev", "N/A"),
        "citation_year": citation.get("year", "N/A"),
        "citation_doi": citation.get("pdbx_database_id_doi", "N/A"),
    }

    return meta


def fetch_polymer_entities(pdb_id: str) -> list:
    """
    Fetch polymer entity details (chains, sequence length, organism).
    """
    pdb_id = pdb_id.upper()
    url = f"{RCSB_API_URL}/polymer_entity/{pdb_id}/1"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    entity_info = {
        "entity_id": 1,
        "description": data.get("rcsb_polymer_entity", {}).get("pdbx_description", "N/A"),
        "organism": data.get("rcsb_entity_source_organism", [{}])[0].get("scientific_name", "N/A")
        if data.get("rcsb_entity_source_organism") else "N/A",
        "sequence_length": data.get("entity_poly", {}).get("rcsb_sample_sequence_length", 0),
        "entity_type": data.get("entity_poly", {}).get("type", "N/A"),
    }

    return [entity_info]


def parse_structure(filepath: Path, pdb_id: str = "PROT"):
    """
    Parse a PDB/mmCIF file into a BioPython Structure object.
    """
    filepath = Path(filepath)
    if filepath.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, str(filepath))


def print_metadata(meta: dict):
    """Pretty-print metadata to console."""
    print("\n" + "=" * 65)
    print(f"  RCSB PDB Entry: {meta['pdb_id']}")
    print("=" * 65)
    print(f"  Title       : {meta['title']}")
    print(f"  Method      : {meta['experimental_method']}")
    print(f"  Resolution  : {meta['resolution']} Å")
    print(f"  Deposited   : {meta['deposition_date'][:10]}")
    print(f"  Mol. Weight : {meta['molecular_weight_kda']} kDa")
    print(f"  Atom Count  : {meta['deposited_atom_count']}")
    print(f"  Residues    : {meta['deposited_residue_count']}")
    print(f"  Disulfides  : {meta['disulfide_bonds']}")
    print(f"  Keywords    : {meta['keywords']}")
    print(f"  Authors     : {', '.join(meta['authors'][:3])}...")
    print(f"  Citation    : {meta['citation_title'][:60]}...")
    print(f"  Journal     : {meta['citation_journal']} ({meta['citation_year']})")
    print(f"  DOI         : {meta['citation_doi']}")
    print("=" * 65 + "\n")
