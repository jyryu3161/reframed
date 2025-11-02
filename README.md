[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI version](https://badge.fury.io/py/reframed.svg)](https://badge.fury.io/py/reframed)
[![Build Status](https://travis-ci.org/cdanielmachado/reframed.svg?branch=master)](https://travis-ci.org/cdanielmachado/reframed)
[![Documentation Status](https://readthedocs.org/projects/reframed/badge/?version=latest)](https://reframed.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/212059108.svg)](https://zenodo.org/badge/latestdoi/212059108)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cdanielmachado/teaching/master?filepath=fba.ipynb)

![ReFramed](reframed_logo.png)

# ReFramed: metabolic modeling package

**ReFramed** implements many constraint-based simulation methods (see list below), and contains interfaces to other
libraries of the *COBRA* ecosystem including [**escher**](https://escher.github.io),
[**cobrapy**](https://opencobra.github.io/cobrapy/), and [**optlang**](https://github.com/biosustain/optlang).

### 💡 New!

A new command-line interface is now available, check the help menu for details:

```
> fba -h
```

### Available methods

Name | Long name | Reference
--- | --- | ---
FBA | Flux Balance Analysis | (Varma and Palsson, 1993)
FVA | Flux Variability Analysis | (Mahadevan and Schilling, 2003)
pFBA | Parsimonious FBA | (Lewis et al, 2010)
FBrAtio | FBA with flux ratios| (Yen et al, 2013)
CAFBA | Constrained Allocation FBA | (Mori et al, 2016)
lMOMA | linear version of MOMA | (Segre et al, 2002)
ROOM | Regulatory On/Off Minimization | (Shlomi et al, 2005)
ll-FBA | loopless FBA | (Schellenberger et al, 2011)
TFA | Thermodynamic Flux Analysis | (Henry et al, 2007)
TVA | Thermodynamic Variability Analysis | (Henry et al, 2007)
NET | Network-embedded thermodynamic analysis |  (Kummel et al 2006)
GIMME | Gene Inactivity Moderated by Met and Exp | (Becker and Palsson, 2008)
E-Flux | E-Flux | (Colijn et al, 2009)
SteadyCom | Community simulation | (Chan et al, 2017)

### 🧬 Transcriptomics Integration (GIMME & E-Flux)

ReFramed provides powerful methods to integrate gene expression data (RNA-seq, microarray) with metabolic models:

#### Quick Start Example

```python
from reframed import load_cbmodel, GIMME, eFlux

# Load model
model = load_cbmodel('ecoli_model.xml')

# Gene expression data from RNA-seq or microarray
gene_expression = {
    'b0008': 10.5,  # Gene ID: expression value
    'b0114': 8.2,
    'b0115': 8.5,
    # ... more genes
}

# GIMME: Penalizes inconsistency with expression data
gimme_solution = GIMME(
    model,
    gene_expression,
    cutoff=25,        # Expression threshold (25th percentile)
    growth_frac=0.9   # Minimum growth (90% of max)
)

# E-Flux: Constrains reaction bounds by expression
eflux_solution = eFlux(model, gene_expression)

print(f"GIMME growth: {gimme_solution.fobj:.4f}")
print(f"E-Flux growth: {eflux_solution.fobj:.4f}")
```

#### Method Comparison

| Feature | GIMME | E-Flux |
|---------|-------|--------|
| **Approach** | Minimizes low-expression reactions | Constrains bounds by expression |
| **Growth constraint** | Explicit minimum growth | Implicit from constraints |
| **Computation** | Two-phase optimization | Single FBA |
| **Robustness** | More robust to noise | Faster, more direct |
| **Best for** | Tissue-specific, cancer metabolism | Large-scale screening |

#### Parameters Guide

**GIMME Parameters:**
- `cutoff` (default: 25): Expression percentile threshold. Lower = more reactions penalized
- `growth_frac` (default: 0.9): Minimum growth fraction. Higher = closer to FBA
- `parsimonious` (default: False): Use parsimonious optimization for simpler flux distribution

**E-Flux Parameters:**
- `scale_rxn` (optional): Reaction ID to normalize fluxes (e.g., glucose uptake)
- `scale_value` (default: 1.0): Scaling factor for normalization
- `parsimonious` (default: False): Use pFBA instead of FBA

#### Complete Examples

See the [`examples/`](examples/) directory for comprehensive examples:

- **`simple_gimme_eflux.py`** - Quick start with E. coli core model
- **`quick_benchmark.py`** - 🆕 Fast benchmark comparing multiple methods
- **`benchmark_transcriptomics.py`** - 🆕 Comprehensive benchmarking framework
- **`custom_method_template.py`** - 🆕 Templates for implementing new methods
- **`ecoli_gimme_eflux_example.py`** - Advanced analysis with:
  - Condition-specific expression simulation (aerobic/anaerobic)
  - Parameter sensitivity analysis
  - Flux comparison between methods
  - Complete workflow demonstration

Run examples:
```bash
cd examples
python simple_gimme_eflux.py          # Basic example
python quick_benchmark.py             # Quick method comparison
python benchmark_transcriptomics.py   # Full benchmark suite
python ecoli_gimme_eflux_example.py   # Detailed analysis
```

#### 🆕 Benchmark Framework

Compare multiple methods systematically:

```python
from benchmark_transcriptomics import TranscriptomicsBenchmark, FBAMethod, GIMMEMethod

benchmark = TranscriptomicsBenchmark(model, gene_expression)
benchmark.add_methods([
    FBAMethod(),
    GIMMEMethod(cutoff=25),
    GIMMEMethod(cutoff=50),
])

results = benchmark.run_all()
benchmark.print_summary()
benchmark.export_results('results.csv')
```

#### 🆕 Adding Custom Methods

Extend the framework with your own algorithms using provided templates:

```python
from benchmark_transcriptomics import TranscriptomicsMethod

class MyMethod(TranscriptomicsMethod):
    def run(self, model, gene_exp, **kwargs):
        # Implement your algorithm
        solution = your_algorithm(model, gene_exp)
        return self._calculate_metrics(solution, self.name, exec_time)

# Use in benchmark
benchmark.add_method(MyMethod())
```

See `custom_method_template.py` for detailed templates and examples.

#### Typical Workflow

1. **Prepare expression data**: Convert to dictionary `{gene_id: expression_value}`
2. **Load metabolic model**: Use `load_cbmodel()` with SBML file
3. **Run analysis**: Apply GIMME or E-Flux
4. **Validate results**: Compare with experimental flux measurements
5. **Perform FVA**: Check flux variability in solution space

```python
from reframed import FVA

# After running GIMME/E-Flux, check flux ranges
variability = FVA(model, obj_frac=0.9)
```

#### References

- **GIMME**: Becker, S. A., & Palsson, B. Ø. (2008). Context-specific metabolic networks are consistent with experiments. *PLoS Computational Biology*, 4(5), e1000082.
- **E-Flux**: Colijn, C., et al. (2009). Interpreting expression data with metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid production. *PLoS Computational Biology*, 5(8), e1000489.

### Documentation

Please check documentation with installation and usage instructions [here](https://reframed.readthedocs.io) or try the live demo [here](https://mybinder.org/v2/gh/cdanielmachado/teaching/master?filepath=fba.ipynb).

Latest version was tested with:

- python 3.13
- libsbml 5.20.4
- PySCIPOpt 5.4.1
- Gurobi 12.0.1
- CPLEX 22.1.1 (with python 3.10)

### Credits and License

Developed by Daniel Machado at the Norwegian University of Science and Technology (2025).

**ReFramed** is a refactored version of the [**framed**](https://github.com/cdanielmachado/framed) library.

Released under an Apache License.

