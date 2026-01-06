# Tracer semi-analytical model (Python)

Python implementation of a modified semi-analytical tracer model used to reproduce and analyze tracer concentration–time behavior using (i) a derived subset of Utah FORGE GDR data and (ii) digitized values from a published curve.

## Data sources

### Utah FORGE (GDR)
The underlying field data are publicly available from the Geothermal Data Repository (GDR):

McLennan, J., England, K., & Swearingen, L. (2024). *Utah FORGE: Wells 16A(78)-32 and 16B(78)-32 Extended Circulation Test Data - August and September 2024* [Data set]. Geothermal Data Repository. https://doi.org/10.15121/2475065

This repository does **not** re-host the full GDR download. The subset used for the analyses is provided as a derived CSV in `data/derived/utah_forge_gdr1683_extracted.csv` and is the default input for the scripts.

### Digitized curve (Rose et al., 2025)
Digitized values from Rose et al. (2025), Figure 1 (NDS concentration vs time), are included in `data/external/rose_2025_fig1_nds_vs_time_digitized.csv` with provenance notes in `data/external/README.md`.

## Repository structure

- `src/` — core model implementation
- `notebooks/` — analysis notebooks (run model, export CSV, generate figures)
- `data/derived/` — derived CSV inputs extracted from the GDR dataset (used by the code)
- `data/external/` — external digitized inputs (published figure data) + provenance
- `data/gdr_raw/` — placeholder for local raw GDR downloads (not required for default runs)
- `docs/` — documentation (data provenance, notes)
- `outputs/` — generated results (CSV outputs). Figures were produced in Excel from these CSV files.

## Requirements

Install required Python packages:

```bash
pip install -r requirements.txt


## Running the notebooks

Install requirements:

```bash
pip install -r requirements.txt

jupyter lab

