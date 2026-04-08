# Protein Structure Visualization Tool

A comprehensive Python toolkit for fetching, analyzing, and visualizing disease-relevant protein structures from the [RCSB Protein Data Bank](https://www.rcsb.org/).

**Target Protein:** SARS-CoV-2 Spike Glycoprotein — [PDB 6VYB](https://www.rcsb.org/structure/6VYB) (open state)

## About the Protein

The **SARS-CoV-2 Spike Glycoprotein (S protein)** is a trimeric class I fusion protein on the surface of the SARS-CoV-2 virus responsible for the COVID-19 pandemic. It mediates viral entry into human cells by binding to the **ACE2 receptor** on host cells.

| Property | Value |
|---|---|
| PDB ID | [6VYB](https://www.rcsb.org/structure/6VYB) |
| Structure | Spike ectodomain (open state) |
| Method | Cryo-EM (3.2 Å resolution) |
| Organism | Severe acute respiratory syndrome coronavirus 2 |
| Molecular Weight | 437.44 kDa |
| Chains | 3 protomers (A, B, C) + glycan chains |
| Atoms | 22,365 deposited |
| Residues | 3,843 (1,281 per protomer) |
| Disulfide Bonds | 35 |
| Citation | Walls et al., *Cell* 181:281 (2020) — [DOI](https://doi.org/10.1016/j.cell.2020.02.058) |

### Key Functional Domains

```
Position 1                                                         1273
|────────────────────────────────────────────────────────────────────|
|  NTD (14-305)  | RBD (319-541) |        |FCS|  FP  | HR1  | CH  |
|                | RBM (437-508) |        |682|      |      |     |
```

- **NTD (N-terminal Domain):** Host cell attachment and immune evasion
- **RBD (Receptor Binding Domain):** Directly contacts the ACE2 receptor; main target for neutralizing antibodies and vaccine design
- **RBM (Receptor Binding Motif):** The specific contact interface with ACE2; contains critical mutation sites (e.g., N501Y, E484K)
- **Furin Cleavage Site (S1/S2):** RRAR motif at position 682-685, unique to SARS-CoV-2 and absent in SARS-CoV — enhances cell entry efficiency
- **Fusion Peptide:** Inserts into the host cell membrane to initiate membrane fusion
- **HR1 / Central Helix:** Structural elements critical for the conformational change that drives virus-cell fusion

## Features

### Static Visualizations (matplotlib)
- **B-factor Profile** — Per-residue thermal displacement with functional region annotations
- **Amino Acid Composition** — Bar chart color-coded by biochemical property
- **Residue Property Distribution** — Donut chart (hydrophobic, polar, charged, special)
- **Cα Distance Matrix (Contact Map)** — Reveals domain boundaries and tertiary contacts
- **Chain B-factor Comparison** — Per-chain flexibility statistics
- **Functional Domain Map** — Linear map of annotated regions

### Interactive 3D Viewers (py3Dmol → HTML)
- **Cartoon View** — Secondary structure colored by spectrum (rainbow N→C)
- **Surface View** — Solvent-accessible surface with cartoon backbone
- **Region Highlights** — Functional domains (RBD, NTD, FP, etc.) in distinct colors
- **B-factor Coloring** — Blue (rigid) to red (flexible) thermal displacement

### Structural Analysis (BioPython)
- Per-chain residue, atom, and B-factor statistics
- Amino acid composition and biochemical property classification
- Cα-Cα distance matrix computation
- Annotated functional regions for known disease proteins
- Full analysis export to JSON

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/protein-viz.git
cd protein-viz
pip install -r requirements.txt
```

### Dependencies

| Package | Purpose |
|---|---|
| `biopython` | PDB parsing, structure analysis |
| `py3Dmol` | Interactive 3D molecular visualization |
| `matplotlib` | Static 2D plots and charts |
| `seaborn` | Enhanced heatmap coloring |
| `numpy` | Numerical computation |
| `requests` | RCSB PDB API and file downloads |
| `Pillow` | Image handling |

## Usage

### Full Pipeline (Default: SARS-CoV-2 Spike)

```bash
python main.py
```

### CLI Options

```bash
python main.py --pdb 6VYB           # Analyze SARS-CoV-2 Spike (default)
python main.py --pdb 6VXX           # Closed-state spike for comparison
python main.py --pdb 7DWZ --chain B # Different protein, different chain
python main.py --info-only           # Print metadata without generating plots
python main.py --no-3d               # Skip interactive HTML generation
```

### As a Python Library

```python
from src.fetcher import fetch_pdb_file, fetch_metadata, parse_structure
from src.analyzer import get_chain_stats, get_residue_composition, summarize_structure
from src.visualizer import plot_bfactor_profile, plot_composition
from src.viewer_3d import create_all_views

# Fetch and parse
pdb_file = fetch_pdb_file("6VYB")
structure = parse_structure(pdb_file, "6VYB")

# Analyze
summary = summarize_structure(structure, "6VYB")
print(f"Total residues: {summary['total_residues']}")
print(f"Total atoms: {summary['total_atoms']}")

# Visualize
metadata = fetch_metadata("6VYB")
from src.analyzer import get_residue_bfactors, identify_key_regions
bfactors = get_residue_bfactors(structure)
regions = identify_key_regions(structure, "6VYB")
plot_bfactor_profile(bfactors, "6VYB", regions)

# Interactive 3D
create_all_views("6VYB")
```

## Output Files

After running the pipeline, the `output/` directory will contain:

```
output/
├── 6VYB_analysis.json      # Full structural analysis (JSON)
├── bfactor_profile.png     # B-factor line plot with region annotations
├── aa_composition.png      # Amino acid bar chart
├── property_dist.png       # Residue property donut chart
├── contact_map.png         # Cα distance matrix heatmap
├── chain_bfactors.png      # Per-chain B-factor comparison
├── domain_map.png          # Linear functional domain map
├── 6VYB_cartoon.html       # 3D cartoon view (open in browser)
├── 6VYB_surface.html       # 3D surface view
├── 6VYB_regions.html       # 3D region highlight view
└── 6VYB_bfactor3d.html     # 3D B-factor colored view
```

Open any `.html` file in a browser for interactive 3D visualization (drag to rotate, scroll to zoom).

## Project Structure

```
protein-viz/
├── main.py                 # CLI entry point
├── requirements.txt        # Python dependencies
├── README.md
├── .gitignore
├── src/
│   ├── __init__.py         # Package metadata
│   ├── fetcher.py          # RCSB PDB download & API queries
│   ├── analyzer.py         # Structural analysis (BioPython)
│   ├── visualizer.py       # Static matplotlib visualizations
│   └── viewer_3d.py        # Interactive 3D HTML viewers (py3Dmol)
├── tests/
│   └── test_pipeline.py    # Unit tests (10 tests)
├── data/                   # Downloaded PDB files (gitignored)
└── output/                 # Generated visualizations (gitignored)
```

## Tests

```bash
python -m unittest tests.test_pipeline -v
```

Runs 10 tests covering:
- PDB file download and caching
- RCSB API metadata fetching
- BioPython structure parsing
- Chain statistics computation
- B-factor extraction
- Amino acid composition analysis
- Biochemical property distribution
- Cα distance matrix computation
- Functional region annotation
- Full structural summary

## Scientific Background

### B-factors (Temperature Factors)

B-factors quantify the displacement of atoms from their mean position in the crystal/cryo-EM structure. Higher B-factors indicate greater atomic mobility or structural disorder. In the spike protein:
- The **RBD** shows elevated B-factors in the open state (6VYB) because it is extended and more flexible
- The **S2 stalk region** (HR1, CH) tends to be more rigid, as it maintains the trimeric architecture

### Contact Maps

The Cα distance matrix reveals the spatial proximity of residue pairs. Off-diagonal clusters of close contacts indicate:
- **Secondary structure** (helices show characteristic diagonal bands)
- **Tertiary contacts** (long-range contacts between distant sequence positions)
- **Domain boundaries** (block-diagonal patterns)

### Why This Protein Matters

PDB 6VYB was among the first cryo-EM structures of the SARS-CoV-2 spike protein, published in *Cell* in February 2020 — just weeks after the virus genome was sequenced. This structure was critical for:
1. Understanding how the virus binds human cells via ACE2
2. Designing mRNA vaccines (Pfizer, Moderna) that encode a stabilized spike
3. Developing monoclonal antibody therapies targeting the RBD
4. Tracking how variants (Alpha, Delta, Omicron) mutate the spike

## Data Source

All structural data is fetched from the [RCSB Protein Data Bank](https://www.rcsb.org/) (RCSB PDB), the primary archive for experimentally determined 3D structures of biological macromolecules.

## License

MIT License

## Author

Naitik Gupta
