# Directionality and Analysis Features

This repository provides an `Analyzer` class and related scripts for directionality studies and general analysis of CYGNO reconstructed tracks.

## Main Components

- **Core Implementation:**  
  - `Analyzer.cxx` (main class implementation)  
  - `Analyzer.h` (class header)

- **Example Scripts:**  
  - `Base_script.cxx`: Basic usage example, suitable for getting started with the analysis workflow.
  - `INAF_dir.cxx`: Advanced implementation supporting command-line parameterization for flexible analysis.
  - `INAF_dir_MP.cxx`: Multiprocessing version compatible with `INAF_dir.cxx`, configurable via a config file for batch processing.

### Mulitiprocessing Usage

The `INAF_dir_MP.cxx` script allows for multiprocessing, it is now the standard for ER directionality analysis.

- **Pile-up Tagging:**  
  The ability to tag pile-up events exists but is currently hardcoded. This feature may be removed or made more flexible in future versions.

- **Event Cuts:**  
  Event cuts are managed through config files; however, the set of variables available for cuts is hardcoded. To introduce new variables, you must update both the config file and the `event_vars` map in the code.

## Python Utilities

- `Analyzer.py`: Minimal Python port of the `Analyzer` class for quick prototyping and analysis.
- `DiffusionAlpha.py`: Example script for alpha diffusion analysis.
- `plotDiffusion.py`: Aggregates results and generates diffusion plots.

## Compilation Instructions

To compile the main analysis program:

```sh
g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum
```

For optimized builds:

```sh
g++ Analyzer.cxx Base_script.cxx -O3 -o nameprog `root-config --libs --cflags` -lSpectrum
```

## Polarimetry Analysis

- `fit_modulation.py`: Script for analyzing results from polarized source measurements.  
  - Fits modulation curves and extracts polarization parameters.

## Analysis Notes & Parameters

- **Rebinning:** Disabled by default in all scripts for raw data fidelity.
- **Recommended Parameters for MANGO He/CF4 Fusion:**
  - `ApplyThr(10)`: Apply threshold of 10 ADC counts.
  - `RemoveNoise(30)`: Remove noise below 30 ADC counts.
  - `NPIP = 250`: Set number of primary ionization points.
  - `wFac = 2.0`: Weighting factor for diffusion.
  - For 17 keV events:  
    `PileUpCandidate(false, counter, true, 0.2, 2.0)`  
    (Disables pile-up, enables candidate search with specified thresholds.)

## Additional Features

- **Track Quality Cuts:**  
  - Scripts include configurable cuts for track length, energy, and cluster size.
- **Output:**  
  - Results are saved in ROOT files for further analysis and plotting.
- **Configurable Parameters:**  
  - Most scripts accept command-line or config file parameters for reproducible workflows.

For further details, refer to comments in each script and the inline documentation in `Analyzer.h`.
