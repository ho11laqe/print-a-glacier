# print-a-glacier

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-research--prototype-orange.svg)

A Python framework for generating 3D-printable STL models of glaciers from NetCDF-based glacier datasets. The tool converts gridded ice thickness and surface elevation fields into physically consistent, watertight meshes suitable for additive manufacturing.

---

## Overview

`print-a-glacier` reconstructs glacier geometry from numerical model output and satellite-derived datasets. It produces two complementary 3D models:

- **Terrain model**: solid topographic base with vertical boundaries and optional elevation offset
- **Glacier model**: watertight ice volume derived from ice thickness and surface elevation fields

The workflow is designed for dual-material 3D printing, typically using opaque PLA for terrain and transparent PLA for glacier ice.

---

## Example Output

![3D printed Aletsch glacier](images/Aletsch3DPrint.jpg)

---

## Scientific Input Data

The tool assumes structured NetCDF input with the following variables:

- `x` — projected x-coordinates  
- `y` — projected y-coordinates  
- `thk` — ice thickness (m)  
- `usurf` — glacier surface elevation (m a.s.l.)

These variables are typically derived from glacier evolution models or OGGM-style simulations.

---

## Methodology

The pipeline consists of the following steps:

1. Load structured NetCDF glacier dataset
2. Compute bedrock topography:  
   \( b = u_{surf} - h_{ice} \)
3. Extract glacier extent using thresholded ice thickness
4. Construct surface triangulation via Delaunay tessellation
5. Enforce watertight geometry by:
   - Filtering long edges
   - Reconstructing boundary walls
6. Export STL geometry for additive manufacturing

---

## Installation

```bash
git clone <repo-url>
cd print-a-glacier
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt