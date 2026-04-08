"""
viewer_3d.py — Interactive 3D protein visualization using py3Dmol.

Generates interactive HTML files with 3D molecular viewers that can be
opened in any web browser. Supports multiple visualization styles:
  - Cartoon (secondary structure colored)
  - Surface (electrostatic-like coloring)
  - Stick (atomic detail)
  - Highlighted functional regions
  - Spike-ACE2 binding interface
"""

import py3Dmol
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent.parent / "output"

# SARS-CoV-2 Spike key regions for highlighting
SPIKE_REGIONS = {
    "RBD": {"start": 319, "end": 541, "color": "#f0883e", "label": "Receptor Binding Domain"},
    "RBM": {"start": 437, "end": 508, "color": "#f85149", "label": "Receptor Binding Motif"},
    "NTD": {"start": 14, "end": 305, "color": "#58a6ff", "label": "N-terminal Domain"},
    "FP": {"start": 788, "end": 806, "color": "#3fb950", "label": "Fusion Peptide"},
    "HR1": {"start": 912, "end": 984, "color": "#d2a8ff", "label": "Heptad Repeat 1"},
    "S1S2": {"start": 682, "end": 685, "color": "#f778ba", "label": "Furin Cleavage Site"},
}


def _wrap_html(viewer_html: str, title: str) -> str:
    """Wrap py3Dmol viewer in a full HTML page with dark theme."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: #0d1117;
            color: #c9d1d9;
            font-family: 'SF Mono', 'Fira Code', monospace;
            display: flex;
            flex-direction: column;
            align-items: center;
            padding: 20px;
            min-height: 100vh;
        }}
        h1 {{
            font-size: 1.4em;
            margin-bottom: 10px;
            color: #58a6ff;
        }}
        .subtitle {{
            font-size: 0.85em;
            color: #8b949e;
            margin-bottom: 20px;
        }}
        .viewer-container {{
            border: 1px solid #30363d;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 4px 24px rgba(0,0,0,0.4);
        }}
        .legend {{
            margin-top: 20px;
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            justify-content: center;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 0.8em;
        }}
        .legend-dot {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
        }}
        .controls {{
            margin-top: 15px;
            font-size: 0.75em;
            color: #8b949e;
        }}
    </style>
</head>
<body>
    <h1>{title}</h1>
    <div class="subtitle">Interactive 3D Viewer — Drag to rotate, scroll to zoom, right-click to translate</div>
    <div class="viewer-container">
        {viewer_html}
    </div>
    <div class="legend" id="legend"></div>
    <div class="controls">
        Built with py3Dmol + BioPython | Data from RCSB PDB
    </div>
</body>
</html>"""


def create_cartoon_view(pdb_id: str = "6VYB", width: int = 900,
                         height: int = 600) -> Path:
    """
    Create an interactive 3D cartoon visualization colored by spectrum.
    Shows secondary structure elements (helices, sheets, loops).
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(f"https://files.rcsb.org/download/{pdb_id}.pdb", "pdb",
                  {"doAssembly": True})

    # Wait - py3Dmol.view() with URL doesn't work in non-notebook mode
    # Instead, we generate the HTML directly with a data fetch
    viewer_html = f"""
    <div id="viewer" style="width: {width}px; height: {height}px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            backgroundColor: "#0d1117"
        }});

        $.ajax({{
            url: "https://files.rcsb.org/download/{pdb_id}.pdb",
            success: function(data) {{
                viewer.addModel(data, "pdb");
                viewer.setStyle({{}}, {{
                    cartoon: {{color: "spectrum", opacity: 0.95}}
                }});
                viewer.zoomTo();
                viewer.spin("y", 0.5);
                viewer.render();
            }}
        }});
    </script>
    """

    html = _wrap_html(viewer_html, f"{pdb_id} — Cartoon View (SARS-CoV-2 Spike)")
    path = OUTPUT_DIR / f"{pdb_id}_cartoon.html"
    path.write_text(html)
    print(f"[3D] Saved cartoon view → {path}")
    return path


def create_surface_view(pdb_id: str = "6VYB", width: int = 900,
                         height: int = 600) -> Path:
    """
    Create a molecular surface view showing the protein's solvent-accessible surface.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    viewer_html = f"""
    <div id="viewer" style="width: {width}px; height: {height}px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            backgroundColor: "#0d1117"
        }});

        $.ajax({{
            url: "https://files.rcsb.org/download/{pdb_id}.pdb",
            success: function(data) {{
                viewer.addModel(data, "pdb");

                // Semi-transparent surface
                viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                    opacity: 0.7,
                    color: "#58a6ff"
                }});

                // Cartoon backbone underneath
                viewer.setStyle({{}}, {{
                    cartoon: {{color: "spectrum", opacity: 0.6}}
                }});

                viewer.zoomTo();
                viewer.spin("y", 0.3);
                viewer.render();
            }}
        }});
    </script>
    """

    html = _wrap_html(viewer_html, f"{pdb_id} — Surface View (SARS-CoV-2 Spike)")
    path = OUTPUT_DIR / f"{pdb_id}_surface.html"
    path.write_text(html)
    print(f"[3D] Saved surface view → {path}")
    return path


def create_region_highlight_view(pdb_id: str = "6VYB", width: int = 900,
                                   height: int = 600) -> Path:
    """
    Highlight functional domains on the spike protein structure.
    RBD, NTD, Fusion Peptide, etc. are shown in distinct colors.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Build style commands for each region
    region_styles = []
    legend_items = []

    for key, region in SPIKE_REGIONS.items():
        region_styles.append(f"""
                viewer.setStyle(
                    {{resi: ["{region['start']}-{region['end']}"], chain: "A"}},
                    {{cartoon: {{color: "{region['color']}", opacity: 1.0, thickness: 0.4}}}}
                );
        """)
        legend_items.append(
            f'<div class="legend-item">'
            f'<div class="legend-dot" style="background:{region["color"]}"></div>'
            f'{region["label"]} ({region["start"]}-{region["end"]})'
            f'</div>'
        )

    styles_js = "\n".join(region_styles)
    legend_html = "\n".join(legend_items)

    viewer_html = f"""
    <div id="viewer" style="width: {width}px; height: {height}px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            backgroundColor: "#0d1117"
        }});

        $.ajax({{
            url: "https://files.rcsb.org/download/{pdb_id}.pdb",
            success: function(data) {{
                viewer.addModel(data, "pdb");

                // Base style: gray cartoon
                viewer.setStyle({{}}, {{
                    cartoon: {{color: "#30363d", opacity: 0.5}}
                }});

                // Highlight regions
                {styles_js}

                viewer.zoomTo();
                viewer.spin("y", 0.4);
                viewer.render();

                // Populate legend
                document.getElementById("legend").innerHTML = `{legend_html}`;
            }}
        }});
    </script>
    """

    html = _wrap_html(viewer_html,
                      f"{pdb_id} — Functional Domain Highlights (SARS-CoV-2 Spike)")
    path = OUTPUT_DIR / f"{pdb_id}_regions.html"
    path.write_text(html)
    print(f"[3D] Saved region highlight view → {path}")
    return path


def create_bfactor_colored_view(pdb_id: str = "6VYB", width: int = 900,
                                  height: int = 600) -> Path:
    """
    Color residues by B-factor (thermal displacement).
    Blue = rigid, Red = flexible.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    viewer_html = f"""
    <div id="viewer" style="width: {width}px; height: {height}px;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            backgroundColor: "#0d1117"
        }});

        $.ajax({{
            url: "https://files.rcsb.org/download/{pdb_id}.pdb",
            success: function(data) {{
                viewer.addModel(data, "pdb");
                viewer.setStyle({{}}, {{
                    cartoon: {{
                        colorscheme: {{
                            prop: "b",
                            gradient: "rwb",
                            min: 0,
                            max: 150
                        }},
                        opacity: 0.95
                    }}
                }});
                viewer.zoomTo();
                viewer.spin("y", 0.4);
                viewer.render();

                document.getElementById("legend").innerHTML = `
                    <div class="legend-item">
                        <div class="legend-dot" style="background:#3333ff"></div>
                        Low B-factor (rigid)
                    </div>
                    <div class="legend-item">
                        <div class="legend-dot" style="background:#ffffff"></div>
                        Medium B-factor
                    </div>
                    <div class="legend-item">
                        <div class="legend-dot" style="background:#ff3333"></div>
                        High B-factor (flexible)
                    </div>
                `;
            }}
        }});
    </script>
    """

    html = _wrap_html(viewer_html, f"{pdb_id} — B-factor Coloring (SARS-CoV-2 Spike)")
    path = OUTPUT_DIR / f"{pdb_id}_bfactor3d.html"
    path.write_text(html)
    print(f"[3D] Saved B-factor 3D view → {path}")
    return path


def create_all_views(pdb_id: str = "6VYB") -> list:
    """Generate all 3D visualization views."""
    views = [
        create_cartoon_view(pdb_id),
        create_surface_view(pdb_id),
        create_region_highlight_view(pdb_id),
        create_bfactor_colored_view(pdb_id),
    ]
    print(f"\n[3D] Generated {len(views)} interactive HTML views")
    return views
