"""
test_pipeline.py — Unit tests for the protein visualization pipeline.
"""

import sys
import json
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.fetcher import fetch_pdb_file, fetch_metadata, parse_structure
from src.analyzer import (
    get_chain_stats, get_residue_bfactors, get_residue_composition,
    get_property_distribution, compute_ca_distance_matrix,
    identify_key_regions, summarize_structure
)


class TestFetcher(unittest.TestCase):
    """Test PDB fetching and parsing."""

    @classmethod
    def setUpClass(cls):
        cls.pdb_id = "6VYB"
        cls.pdb_file = fetch_pdb_file(cls.pdb_id)
        cls.structure = parse_structure(cls.pdb_file, cls.pdb_id)

    def test_file_downloaded(self):
        self.assertTrue(self.pdb_file.exists())
        self.assertGreater(self.pdb_file.stat().st_size, 1000)

    def test_metadata(self):
        meta = fetch_metadata(self.pdb_id)
        self.assertEqual(meta["pdb_id"], "6VYB")
        self.assertIn("spike", meta["title"].lower())
        self.assertIsInstance(meta["authors"], list)
        self.assertGreater(len(meta["authors"]), 0)
        self.assertGreater(meta["deposited_atom_count"], 0)

    def test_structure_parsed(self):
        model = self.structure[0]
        chains = list(model.get_chains())
        self.assertGreater(len(chains), 0)


class TestAnalyzer(unittest.TestCase):
    """Test structural analysis functions."""

    @classmethod
    def setUpClass(cls):
        cls.pdb_id = "6VYB"
        pdb_file = fetch_pdb_file(cls.pdb_id)
        cls.structure = parse_structure(pdb_file, cls.pdb_id)

    def test_chain_stats(self):
        stats = get_chain_stats(self.structure)
        self.assertGreater(len(stats), 0)
        for cs in stats:
            self.assertIn("chain_id", cs)
            self.assertIn("residue_count", cs)
            self.assertGreater(cs["atom_count"], 0)

    def test_bfactors(self):
        bf = get_residue_bfactors(self.structure)
        self.assertGreater(len(bf["bfactors"]), 0)
        self.assertEqual(len(bf["bfactors"]), len(bf["residue_ids"]))
        self.assertEqual(len(bf["bfactors"]), len(bf["residue_names"]))

    def test_composition(self):
        comp = get_residue_composition(self.structure)
        self.assertGreater(comp["total_residues"], 0)
        self.assertIn("composition", comp)
        # Should have standard amino acids
        self.assertGreater(len(comp["composition"]), 10)

    def test_property_distribution(self):
        props = get_property_distribution(self.structure)
        self.assertIn("distribution", props)
        self.assertIn("Hydrophobic", props["distribution"])
        total_pct = sum(
            v["percentage"] for v in props["distribution"].values()
        )
        self.assertAlmostEqual(total_pct, 100.0, places=0)

    def test_distance_matrix(self):
        dist = compute_ca_distance_matrix(self.structure, "A", max_residues=50)
        self.assertEqual(dist["distance_matrix"].shape[0], dist["n_residues"])
        self.assertEqual(dist["distance_matrix"].shape[1], dist["n_residues"])
        # Diagonal should be zero
        for i in range(dist["n_residues"]):
            self.assertAlmostEqual(dist["distance_matrix"][i, i], 0.0, places=3)

    def test_key_regions(self):
        regions = identify_key_regions(self.structure, "6VYB")
        self.assertGreater(len(regions), 0)
        for r in regions:
            self.assertIn("name", r)
            self.assertIn("start", r)
            self.assertIn("end", r)
            self.assertGreater(r["end"], r["start"])

    def test_summarize(self):
        summary = summarize_structure(self.structure, "6VYB")
        self.assertGreater(summary["total_atoms"], 0)
        self.assertGreater(summary["total_residues"], 0)
        self.assertGreater(summary["total_chains"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
