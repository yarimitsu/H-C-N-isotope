# Bayesian Three-Isotope Mixing Model for Glacier-Marine Food Webs

This repository contains the code and data for a novel Bayesian three-isotope (C-N-H) mixing model that uses deuterium (δ²H) as a third isotope alongside carbon and nitrogen to trace biogeochemical subsidies from glacier runoff into Alaska's coastal marine food webs.

## Publication

**Arimitsu, M. L., Hobson, K. A., Webber, D. N., Piatt, J. F., Hood, E. W., & Fellman, J. B. (2018).** Tracing biogeochemical subsidies from glacier runoff into Alaska's coastal marine food webs. *Global Change Biology*, 24(1), 387-398. https://doi.org/10.1111/gcb.13875

**USGS Publication:** https://www.usgs.gov/publications/tracing-biogeochemical-subsidies-glacier-runoff-alaskas-coastal-marine-food-webs

## Overview

This project implements a multi-trophic level Bayesian three-isotope mixing model that incorporates deuterium (δ²H) alongside traditional carbon (δ¹³C) and nitrogen (δ¹⁵N) stable isotopes**, plus radiogenic carbon isotopes (Δ¹⁴C), to trace riverine and marine organic matter sources as they are transferred through coastal marine food webs in glacier-influenced habitats of the Gulf of Alaska.

### Methodological Innovations

#### 1. Deuterium as a Third Isotope

The inclusion of hydrogen isotopes (δ²H) represents a key methodological advancement in food web ecology. While most isotope mixing models rely on carbon and nitrogen isotopes alone, this study leverages the large natural gradients in δ²H between glacial freshwater sources and marine waters. This third isotope dimension:

- **Provides stronger source discrimination** between marine, coastal, and freshwater organic matter inputs
- **Accounts for the dual origin of tissue hydrogen**: both dietary sources and ambient water (via drinking water and metabolic water production)
- **Enables more robust estimates** of freshwater subsidies to marine food webs than traditional two-isotope approaches
- **Exploits the unique isotopic signature** of glacier-fed rivers, which are substantially depleted in deuterium compared to seawater

The model explicitly partitions hydrogen isotope ratios in consumer tissues into contributions from diet and ambient water using the parameter ω (omega), which varies by species based on trophic level and physiology.

#### 2. Taxa-Tissue Dependent Trophic Discrimination

The model incorporates tissue-specific trophic discrimination factors, recognizing that isotopic fractionation varies among tissue types. Rather than applying a single fractionation factor across all organisms, the model estimates separate discrimination factors for six tissue/taxa categories:

- **Bird liver** (Δ¹³C, Δ¹⁵N)
- **Bird blood** (Δ¹³C, Δ¹⁵N)
- **Fish muscle** (Δ¹³C, Δ¹⁵N)
- **Macrozooplankton** (Δ¹³C, Δ¹⁵N)
- **Microzooplankton** (Δ¹³C, Δ¹⁵N)
- **Young-of-year fish** (Δ¹³C, Δ¹⁵N)

This approach accounts for the well-documented biological reality that trophic fractionation differs between tissue types (e.g., liver vs. blood in birds, mussel vs. whole body tissue in fish) and among taxonomic groups (birds, fish, zooplankton). By incorporating these taxa and tissue-specific parameters, the model produces more accurate estimates of trophic position and source contributions than models using uniform fractionation factors.

### Key Findings

- Riverine organic matter assimilation ranged from 12% to 44% across species at higher trophic levels
- Different species showed varying reliance on terrestrial subsidies based on their habitat preferences and life histories
- While organic carbon in glacial runoff was ancient (12,100–1,500 years old), coastal organisms utilized much younger material (530 years to modern age)
- This suggests ecosystems are processing recently mobilized terrestrial carbon rather than deeply aged deposits

## Model Description

The model uses a Bayesian framework implemented in Stan to estimate:

1. **Source contributions** (offshore marine, coastal, and freshwater) to consumer diets across three isotope systems simultaneously (δ¹³C, δ¹⁵N, δ²H)
2. **Trophic fractionation** for carbon, nitrogen, and hydrogen isotopes across different tissue types
3. **Trophic levels** for each species group
4. **Proportion of tissue hydrogen (ω)** derived from ambient water vs. diet — a critical parameter unique to hydrogen isotope modeling

### Novel Features for Hydrogen Isotopes

Unlike carbon and nitrogen, hydrogen in animal tissues has two sources:
- **Dietary organic hydrogen** (from consumed prey)
- **Ambient water hydrogen** (from drinking water and metabolic water production)

The model addresses this by:
- Estimating species-specific ω (omega) parameters that partition tissue δ²H between these sources
- Calculating diet-derived δ²H values (δ²H_diet) by correcting for ambient water contributions
- Accounting for trophic-level dependent changes in ω for higher consumers
- Using measured δ²H of ambient water as a covariate in the model

The model also accounts for:
- Multiple species groups (17 species/taxa)
- Different tissue types (bird liver, bird blood, fish, macrozooplankton, microzooplankton, young-of-year fish)
- Three potential organic matter sources (offshore, coastal, freshwater)
- Isotopic fractionation at each trophic level for all three isotope systems

## Repository Structure

```
H-C-N-isotope/
├── Main.R              # Main script to run the Bayesian MCMC model
├── Plot.R              # Script to generate diagnostic and results plots
├── isotope7.stan       # Stan model specification
├── data/
│   └── glacier_iso.csv # Input isotope data
├── figs/               # Output directory for generated plots
└── README.md           # This file
```

## Data

The input data (`data/glacier_iso.csv`) contains isotope measurements for various marine organisms collected from Gulf of Alaska glacier-marine habitats, including:

- **Taxa**: Species taxonomic group and tissue type
- **Species**: Species identifier
- **Isotope measurements**: δ¹³C, δ¹⁵N, δ²H
- **Environmental covariates**: Year, fjord location, water δ²H
- **Trophic level estimates**

Species analyzed include:
- Seabirds: Black-legged Kittiwakes (BLKI), Kittlitz's Murrelets (KIMU), Marbled Murrelets (MAMU)
- Forage fish: Capelin, Herring, Eulachon, Sandlance, Age-0 Walleye Pollock 
- Zooplankton: Copepods, Euphausiids (*E. pacifica*, *Thysanoessa*), Amphipods (*Themisto libellula*), Mysids (*Neomysis rayii*)

## Requirements

The code requires R with the following packages:

```r
rstan          # Bayesian inference using Stan
reshape2       # Data manipulation
ggplot2        # Plotting
plyr           # Data manipulation
png            # PNG graphics device
scales         # Plotting utilities
```

## Usage

1. **Ensure all dependencies are installed:**
   ```r
   install.packages(c("rstan", "reshape2", "ggplot2", "plyr", "png", "scales"))
   ```

2. **Run the main analysis:**
   ```r
   source("Main.R")
   ```

   This script will:
   - Load the isotope data
   - Compile the Stan model
   - Run MCMC sampling (100,000 iterations with 10,000 warmup, thinning by 90)
   - Extract and summarize posterior distributions
   - Generate diagnostic plots via `Plot.R`

3. **Computational requirements:**
   - The MCMC sampling is computationally intensive and may take several hours to complete
   - Multi-core processing is enabled by default using `parallel::detectCores()`

## Output

The analysis generates:

1. **MCMC diagnostics:**
   - Trace plots for model parameters
   - Posterior distributions with prior overlays
   - Convergence diagnostics

2. **Results plots** (saved in `figs/`):
   - Isotope biplots (δ¹³C vs δ²H, δ¹³C vs δ¹⁵N, δ¹⁵N vs δ²H)
   - Source contribution posteriors by species
   - Trophic level estimates
   - Model fit visualization with posterior predictive distributions
   - Fractionation corrected isotope mixing spaces

3. **Summary statistics:**
   - Posterior medians and standard deviations for all parameters
   - Source proportion estimates (phi) for each species
   - Trophic level estimates (tau) for each species

## Model Parameters

Key parameters estimated by the model:

**Source and Diet Composition:**
- `phi[j,k]`: Proportion of source k (offshore/coastal/freshwater) in the diet of species j (simplex constraint: sums to 1)
- `dC1, dC2, dC3`: δ¹³C values of offshore, coastal, and freshwater sources
- `dN[k]`, `dH[k]`: δ¹⁵N and δ²H values for each source (k = 1:3)

**Trophic Structure:**
- `tau[j]`: Trophic level of species j (constrained between 1.8 and 5.5)
- `dN_base`: δ¹⁵N value of primary consumers (baseline for trophic level calculations)

**Taxa-Tissue Specific Trophic Fractionation:**
- `Delta_C[m]`: Carbon trophic discrimination factor for tissue type m (m = 1:6)
- `Delta_N[m]`: Nitrogen trophic discrimination factor for tissue type m (m = 1:6)
  - m=1: bird liver, m=2: bird blood, m=3: fish, m=4: macrozooplankton, m=5: microzooplankton, m=6: YOY fish
- `Delta_H`: Hydrogen trophic fractionation (water-to-phytoplankton)

**Hydrogen Isotope Partitioning (Novel):**
- `omega[j]`: Proportion of tissue hydrogen derived from ambient water (species-specific)
- `dH_w[j]`: Measured δ²H of ambient water (covariate)
- `fblki`: Additional hydrogen fractionation factor for Black-legged Kittiwakes

**Uncertainty Parameters:**
- `sigma_src[i,k]`: Standard deviation of source isotope values (i = C, N, H; k = offshore, coastal, freshwater)
- `sigma_frc[i,j]`: Standard deviation of fractionation by isotope system and species
- `sigma_frcH`, `sigma_omega`, `sigma_fblki`: Additional variance parameters for hydrogen modeling

## Citation

If you use this code or data, please cite:

Arimitsu, M. L., Hobson, K. A., Webber, D. N., Piatt, J. F., Hood, E. W., & Fellman, J. B. (2018). Tracing biogeochemical subsidies from glacier runoff into Alaska's coastal marine food webs. *Global Change Biology*, 24(1), 387-398. https://doi.org/10.1111/gcb.13875

## Contact

For questions about the methodology or data, please contact:
- Mayumi Arimitsu, Ocean Bight LLC, marimitsu@oceanbight.com

## License

This project is maintained as a scientific research archive.

## Acknowledgments

This research was funded by the North Pacific Research Board and USGS Alaska Science Center to examine how climate-driven changes in freshwater discharge could reshape Alaska's coastal food webs.
