# pyexpdes

Python implementation of examples from *Experimental Design for Data Science and Engineering* by V. Roshan Joseph.

## Installation

```bash
pip install -e .
```

## Project Structure

```
pyexpdes/
├── src/pyexpdes/
│   ├── core/           # paths, style, R wrapper utilities
│   ├── chapters/       # ch01-ch11 modules
│   └── r_wrappers/     # R scripts by chapter (ch01-ch11)
├── notebooks/          # Jupyter notebooks (ch01-ch11)
├── datasets/           # Data files
└── figures/            # Generated figures (ch01-ch11)
```

## Usage

```python
from pyexpdes.chapters.ch01_introduction import generate_figures

out_dir, fig_paths = generate_figures(use_r=True)
```

Or run notebooks:

```bash
jupyter lab notebooks/ch01.ipynb
```

## Chapters

| Module | Topic |
|--------|-------|
| ch01_introduction | Introduction |
| ch02_modeling | Modeling |
| ch03_model_based_designs | Model-based Designs |
| ch04_space_filling_designs | Space-Filling Designs |
| ch05_representative_points | Representative Points |
| ch06_screening_designs | Screening Designs |
| ch07_sequential_designs | Sequential Designs |
| ch08_multilevel_mixture | Fractional Factorial Designs |
| ch09_model_calibration | Model Calibration |
| ch10_data_subsampling | Data Subsampling |
| ch11_data_analysis | Data Analysis |

## Dependencies

### Python
- numpy, scipy, pandas, matplotlib

### R (for wrapper-backed computations)
- Rscript must be on PATH

**Core:** support, rkriging, AlgDesign

**By chapter:**
- Ch4: minimaxdesign, MaxPro, SFDesign, mined, maximin, SLHD, fitdistrplus, cubature, adaptMCMC
- Ch5: SFDesign, spacefillr, mnormt, adaptMCMC
- Ch6: sensitivity, MOFAT
- Ch7: MaxPro, OSFD, FNN, pdist
- Ch8: HiGarrote, DoE.base, twinning, mixexp
- Ch9: ICAOD, numDeriv
- Ch10: SPlit, BalancedSampling, supercompress
- Ch11: first, randomForest, twingp

`support`, `mixexp`, `minimaxdesign` and `mined` have been archived on CRAN,
so `install.packages()` no longer finds them. They install from the archive:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/support/support_0.1.7.tar.gz",
  repos = NULL, type = "source"
)
```

The other three are under `Archive/mixexp/`, `Archive/minimaxdesign/` and
`Archive/mined/` at the same address. Building them needs a C++ / Fortran
toolchain.

## Related

- R implementation: [vroshanjoseph/rexpdes](https://github.com/vroshanjoseph/rexpdes)

## Notes

- zipcodeR excluded due to package conflicts
- Ch5 PQMC section (Figure 5.9) excluded due to external dependency
