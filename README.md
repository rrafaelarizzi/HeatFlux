# Heat Flux Calculation and Submesoscale Filtering

This repository contains two MATLAB routines developed to compute and isolate latent and sensible heat fluxes associated with submesoscale oceanic processes.

## 🧮 1. `calculate_LHF_SHF.m` — Latent and Sensible Heat Flux Calculation

This script calculates **latent heat flux (LHF)** and **sensible heat flux (SHF)** using daily atmospheric and oceanographic data.

### 🔧 Input Variables
- Air temperature (°C)
- Dew point temperature (°C)
- Wind speed (m/s)
- Air density (kg/m³)
- Sea surface temperature (°C)

### ⚙️ Output
- Latent Heat Flux (LHF) in W/m²
- Sensible Heat Flux (SHF) in W/m²

These fluxes represent the total vertical heat transfer at the ocean–atmosphere interface.

## 🌀 2. `Filtering_Fluxes_Hann.m` — Submesoscale-Associated Flux Extraction

This script applies **Hann filters** in both time and space to isolate the portion of the LHF and SHF associated with **submesoscale features**, such as oceanic vortices.

### 🔧 Filtering Procedure
- **Temporal filtering**: isolates processes with periods shorter than a chosen cutoff.
- **Spatial filtering**: isolates horizontal scales typically associated with submesoscale variability.

### ⚙️ Output
- Filtered LHF and SHF fields representing fluxes linked to submesoscale activity.

## 📁 File Structure

```
├── calculate_LHF_SHF.m          # Computes LHF and SHF from daily data
├── Filtering_Fluxes_Hann.m      # Applies filters to isolate submesoscale-associated fluxes
├── README.md                    # This file
```

## 👩‍💻 Author

Rafaela Rizzi — [rafaela.rizzi@gmail.com](mailto:rafaela.rizzi@gmail.com)

2022–2025

---

Feel free to cite or reuse with credit. Contributions and suggestions are welcome!
