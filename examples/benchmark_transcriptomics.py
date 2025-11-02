#!/usr/bin/env python
"""
Transcriptomics Integration Methods Benchmark

This script provides a comprehensive benchmarking framework for comparing
different methods of integrating gene expression data with metabolic models.

Features:
- Compare multiple methods with same data
- Automatic performance metrics calculation
- Results visualization and export
- Easy extension with custom methods

Author: ReFramed
Date: 2025
"""

import time
import os
import numpy as np
import pandas as pd
from typing import Dict, List, Callable, Optional, Tuple, Any
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from reframed import load_cbmodel, FBA, pFBA, GIMME, eFlux
from reframed.solvers.solution import Status


@dataclass
class BenchmarkResult:
    """Store results from a single method execution."""
    method_name: str
    status: str
    growth_rate: float
    execution_time: float
    active_reactions: int
    total_flux: float
    flux_values: Dict[str, float] = field(default_factory=dict)
    parameters: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None


class TranscriptomicsMethod(ABC):
    """
    Abstract base class for transcriptomics integration methods.

    Extend this class to add new methods to the benchmark.
    """

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        """
        Execute the method and return results.

        Parameters
        ----------
        model : CBModel
            The metabolic model
        gene_exp : dict
            Gene expression data {gene_id: expression_value}
        **kwargs
            Additional method-specific parameters

        Returns
        -------
        BenchmarkResult
            Results including growth rate, fluxes, and metrics
        """
        pass

    @staticmethod
    def _calculate_metrics(solution, method_name: str, exec_time: float,
                          parameters: Dict = None) -> BenchmarkResult:
        """Helper to calculate standard metrics from solution."""
        if parameters is None:
            parameters = {}

        if solution.status.name != 'OPTIMAL':
            return BenchmarkResult(
                method_name=method_name,
                status=solution.status.name,
                growth_rate=0.0,
                execution_time=exec_time,
                active_reactions=0,
                total_flux=0.0,
                flux_values={},
                parameters=parameters,
                error=f"Non-optimal status: {solution.status.name}"
            )

        flux_values = solution.values
        active_rxns = sum(1 for v in flux_values.values() if abs(v) > 1e-6)
        total_flux = sum(abs(v) for v in flux_values.values())

        return BenchmarkResult(
            method_name=method_name,
            status=solution.status.name,
            growth_rate=solution.fobj,
            execution_time=exec_time,
            active_reactions=active_rxns,
            total_flux=total_flux,
            flux_values=dict(flux_values),
            parameters=parameters
        )


class FBAMethod(TranscriptomicsMethod):
    """Standard FBA (baseline - no expression data used)."""

    def __init__(self):
        super().__init__("FBA")

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        start_time = time.time()
        solution = FBA(model)
        exec_time = time.time() - start_time

        return self._calculate_metrics(solution, self.name, exec_time)


class pFBAMethod(TranscriptomicsMethod):
    """Parsimonious FBA (baseline - no expression data used)."""

    def __init__(self):
        super().__init__("pFBA")

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        start_time = time.time()
        solution = pFBA(model)
        exec_time = time.time() - start_time

        return self._calculate_metrics(solution, self.name, exec_time)


class GIMMEMethod(TranscriptomicsMethod):
    """GIMME method with configurable parameters."""

    def __init__(self, cutoff: float = 25, growth_frac: float = 0.9,
                 parsimonious: bool = False):
        super().__init__(f"GIMME(cutoff={cutoff},gf={growth_frac})")
        self.cutoff = cutoff
        self.growth_frac = growth_frac
        self.parsimonious = parsimonious

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        try:
            start_time = time.time()
            solution = GIMME(
                model,
                gene_exp,
                cutoff=self.cutoff,
                growth_frac=self.growth_frac,
                parsimonious=self.parsimonious
            )
            exec_time = time.time() - start_time

            params = {
                'cutoff': self.cutoff,
                'growth_frac': self.growth_frac,
                'parsimonious': self.parsimonious
            }

            return self._calculate_metrics(solution, self.name, exec_time, params)
        except Exception as e:
            return BenchmarkResult(
                method_name=self.name,
                status="ERROR",
                growth_rate=0.0,
                execution_time=0.0,
                active_reactions=0,
                total_flux=0.0,
                parameters={'cutoff': self.cutoff, 'growth_frac': self.growth_frac},
                error=str(e)
            )


class EFluxMethod(TranscriptomicsMethod):
    """E-Flux method with configurable parameters."""

    def __init__(self, scale_rxn: Optional[str] = None, scale_value: float = 1.0,
                 parsimonious: bool = False):
        name_parts = ["E-Flux"]
        if scale_rxn:
            name_parts.append(f"scale={scale_rxn}")
        super().__init__("_".join(name_parts))

        self.scale_rxn = scale_rxn
        self.scale_value = scale_value
        self.parsimonious = parsimonious

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        try:
            start_time = time.time()
            solution = eFlux(
                model,
                gene_exp,
                scale_rxn=self.scale_rxn,
                scale_value=self.scale_value,
                parsimonious=self.parsimonious
            )
            exec_time = time.time() - start_time

            params = {
                'scale_rxn': self.scale_rxn,
                'scale_value': self.scale_value,
                'parsimonious': self.parsimonious
            }

            return self._calculate_metrics(solution, self.name, exec_time, params)
        except Exception as e:
            return BenchmarkResult(
                method_name=self.name,
                status="ERROR",
                growth_rate=0.0,
                execution_time=0.0,
                active_reactions=0,
                total_flux=0.0,
                parameters={'scale_rxn': self.scale_rxn},
                error=str(e)
            )


class TranscriptomicsBenchmark:
    """
    Benchmark framework for comparing transcriptomics integration methods.

    Example usage:
    >>> benchmark = TranscriptomicsBenchmark(model, gene_expression)
    >>> benchmark.add_method(FBAMethod())
    >>> benchmark.add_method(GIMMEMethod(cutoff=25, growth_frac=0.9))
    >>> results = benchmark.run_all()
    >>> benchmark.print_summary()

    With essentiality evaluation:
    >>> benchmark.set_essentiality_data(essentiality_data)
    >>> benchmark.run_all_with_essentiality()
    >>> benchmark.print_essentiality_summary()
    """

    def __init__(self, model, gene_expression: Dict[str, float],
                 essentiality_data: Optional[Dict] = None):
        """
        Initialize benchmark with model and expression data.

        Parameters
        ----------
        model : CBModel
            The metabolic model to benchmark
        gene_expression : dict
            Gene expression data {gene_id: expression_value}
        essentiality_data : dict, optional
            Gene essentiality data for model validation
        """
        self.model = model
        self.gene_expression = gene_expression
        self.methods: List[TranscriptomicsMethod] = []
        self.results: List[BenchmarkResult] = []
        self.essentiality_data = essentiality_data
        self.essentiality_evaluations: Dict[str, Any] = {}

    def add_method(self, method: TranscriptomicsMethod):
        """Add a method to the benchmark."""
        self.methods.append(method)

    def add_methods(self, methods: List[TranscriptomicsMethod]):
        """Add multiple methods to the benchmark."""
        self.methods.extend(methods)

    def run_all(self, verbose: bool = True) -> List[BenchmarkResult]:
        """
        Run all added methods and collect results.

        Parameters
        ----------
        verbose : bool
            Print progress messages

        Returns
        -------
        list of BenchmarkResult
            Results from all methods
        """
        self.results = []

        if verbose:
            print(f"\n{'='*70}")
            print(f"Running Benchmark: {len(self.methods)} methods")
            print(f"Model: {len(self.model.reactions)} reactions, "
                  f"{len(self.model.metabolites)} metabolites, "
                  f"{len(self.model.genes)} genes")
            print(f"Expression data: {len(self.gene_expression)} genes")
            print(f"{'='*70}\n")

        for i, method in enumerate(self.methods, 1):
            if verbose:
                print(f"[{i}/{len(self.methods)}] Running {method.name}...", end=" ")

            try:
                result = method.run(self.model, self.gene_expression)
                self.results.append(result)

                if verbose:
                    if result.error:
                        print(f"✗ ERROR: {result.error}")
                    else:
                        print(f"✓ Growth: {result.growth_rate:.4f}, "
                              f"Time: {result.execution_time:.3f}s")
            except Exception as e:
                error_result = BenchmarkResult(
                    method_name=method.name,
                    status="EXCEPTION",
                    growth_rate=0.0,
                    execution_time=0.0,
                    active_reactions=0,
                    total_flux=0.0,
                    error=str(e)
                )
                self.results.append(error_result)

                if verbose:
                    print(f"✗ EXCEPTION: {str(e)}")

        return self.results

    def print_summary(self):
        """Print a formatted summary of benchmark results."""
        if not self.results:
            print("No results available. Run benchmark first.")
            return

        print(f"\n{'='*100}")
        print("BENCHMARK SUMMARY")
        print(f"{'='*100}")

        # Header
        print(f"{'Method':<35} {'Status':<10} {'Growth':>10} {'Active Rxns':>12} "
              f"{'Total Flux':>12} {'Time (s)':>10}")
        print("-" * 100)

        # Results
        for result in self.results:
            status_symbol = "✓" if result.status == "OPTIMAL" else "✗"
            print(f"{result.method_name:<35} {status_symbol} {result.status:<8} "
                  f"{result.growth_rate:>10.4f} {result.active_reactions:>12} "
                  f"{result.total_flux:>12.2f} {result.execution_time:>10.3f}")

        print("=" * 100)

    def compare_fluxes(self, reaction_ids: List[str], top_n: int = None):
        """
        Compare flux values for specific reactions across methods.

        Parameters
        ----------
        reaction_ids : list of str
            Reaction IDs to compare
        top_n : int, optional
            If provided, show only top N reactions by variance
        """
        if not self.results:
            print("No results available. Run benchmark first.")
            return

        print(f"\n{'='*100}")
        print("FLUX COMPARISON")
        print(f"{'='*100}")

        # Create DataFrame for easy comparison
        flux_data = {}
        for result in self.results:
            if result.status == "OPTIMAL":
                flux_data[result.method_name] = result.flux_values

        if not flux_data:
            print("No optimal solutions to compare.")
            return

        # Determine reactions to show
        if top_n:
            # Calculate variance across methods for each reaction
            all_rxns = set()
            for fluxes in flux_data.values():
                all_rxns.update(fluxes.keys())

            rxn_variance = {}
            for rxn in all_rxns:
                values = [flux_data[method].get(rxn, 0)
                         for method in flux_data.keys()]
                if any(abs(v) > 1e-6 for v in values):  # Only consider active reactions
                    rxn_variance[rxn] = np.var(values)

            # Get top N by variance
            sorted_rxns = sorted(rxn_variance.items(), key=lambda x: x[1], reverse=True)
            reaction_ids = [rxn for rxn, _ in sorted_rxns[:top_n]]

            print(f"Showing top {top_n} reactions by variance")

        # Print comparison
        method_names = list(flux_data.keys())
        header = f"{'Reaction':<20}"
        for name in method_names:
            header += f"{name[:18]:>20}"
        print(header)
        print("-" * (20 + 20 * len(method_names)))

        for rxn_id in reaction_ids:
            row = f"{rxn_id:<20}"
            for method_name in method_names:
                flux_val = flux_data[method_name].get(rxn_id, 0)
                row += f"{flux_val:>20.4f}"
            print(row)

        print("=" * 100)

    def export_results(self, filename: str, format: str = 'csv'):
        """
        Export benchmark results to file.

        Parameters
        ----------
        filename : str
            Output filename
        format : str
            Output format: 'csv', 'json', or 'excel'
        """
        if not self.results:
            print("No results available. Run benchmark first.")
            return

        # Create summary DataFrame
        data = []
        for result in self.results:
            row = {
                'Method': result.method_name,
                'Status': result.status,
                'Growth_Rate': result.growth_rate,
                'Active_Reactions': result.active_reactions,
                'Total_Flux': result.total_flux,
                'Execution_Time': result.execution_time,
                'Error': result.error or ''
            }
            # Add parameters
            for key, value in result.parameters.items():
                row[f'Param_{key}'] = value
            data.append(row)

        df = pd.DataFrame(data)

        if format == 'csv':
            df.to_csv(filename, index=False)
        elif format == 'json':
            df.to_json(filename, orient='records', indent=2)
        elif format == 'excel':
            df.to_excel(filename, index=False)
        else:
            raise ValueError(f"Unknown format: {format}")

        print(f"\nResults exported to: {filename}")

    def export_detailed_fluxes(self, filename: str):
        """
        Export detailed flux values for all methods and reactions.

        Parameters
        ----------
        filename : str
            Output CSV filename
        """
        if not self.results:
            print("No results available. Run benchmark first.")
            return

        # Collect all reactions
        all_reactions = set()
        for result in self.results:
            all_reactions.update(result.flux_values.keys())

        # Create flux comparison table
        data = {'Reaction': sorted(all_reactions)}
        for result in self.results:
            data[result.method_name] = [
                result.flux_values.get(rxn, 0.0) for rxn in data['Reaction']
            ]

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)

        print(f"\nDetailed fluxes exported to: {filename}")

    def set_essentiality_data(self, essentiality_data: Dict):
        """
        Set gene essentiality data for model validation.

        Parameters
        ----------
        essentiality_data : dict
            Gene essentiality data from load_essentiality_data()
        """
        self.essentiality_data = essentiality_data
        print(f"✓ Essentiality data set for {len(essentiality_data)} genes")

    def evaluate_with_essentiality(self, growth_threshold: float = 0.01) -> pd.DataFrame:
        """
        Evaluate all methods using gene essentiality data.

        Parameters
        ----------
        growth_threshold : float
            Growth rate threshold for essentiality prediction

        Returns
        -------
        DataFrame
            Essentiality evaluation results
        """
        if self.essentiality_data is None:
            raise ValueError("Essentiality data not set. Use set_essentiality_data() first.")

        # Import here to avoid circular dependency
        try:
            from gene_essentiality_evaluation import GeneEssentialityEvaluator
        except ImportError:
            print("Error: gene_essentiality_evaluation module not found")
            return None

        evaluator = GeneEssentialityEvaluator(
            self.model,
            self.essentiality_data,
            growth_threshold=growth_threshold
        )

        eval_results = []

        print(f"\n{'='*70}")
        print("Evaluating Methods with Gene Essentiality Data")
        print(f"{'='*70}")

        for result in self.results:
            method_name = result.method_name

            # Determine method type
            if 'FBA' in method_name and 'pFBA' not in method_name and 'GIMME' not in method_name:
                predictions = evaluator.predict_essentiality(method='fba')
            elif 'pFBA' in method_name:
                predictions = evaluator.predict_essentiality(method='pfba')
            elif 'GIMME' in method_name:
                # Extract parameters from method_name if possible
                params = result.parameters
                predictions = evaluator.predict_essentiality(
                    method='gimme',
                    gene_exp=self.gene_expression,
                    **params
                )
            elif 'E-Flux' in method_name or 'eFlux' in method_name:
                params = result.parameters
                predictions = evaluator.predict_essentiality(
                    method='eflux',
                    gene_exp=self.gene_expression,
                    **params
                )
            else:
                print(f"Skipping {method_name}: Unknown method type")
                continue

            evaluation = evaluator.calculate_metrics(predictions, method_name)
            self.essentiality_evaluations[method_name] = evaluation

            eval_results.append({
                'Method': method_name,
                'Accuracy': evaluation.accuracy,
                'Precision': evaluation.precision,
                'Recall': evaluation.recall,
                'F1-Score': evaluation.f1_score,
                'MCC': evaluation.mcc,
                'TP': evaluation.true_positives,
                'TN': evaluation.true_negatives,
                'FP': evaluation.false_positives,
                'FN': evaluation.false_negatives,
            })

            print(f"✓ {method_name}: Accuracy={evaluation.accuracy:.3f}, "
                  f"F1={evaluation.f1_score:.3f}")

        return pd.DataFrame(eval_results)

    def print_essentiality_summary(self):
        """Print summary of essentiality evaluation results."""
        if not self.essentiality_evaluations:
            print("No essentiality evaluations available.")
            print("Run evaluate_with_essentiality() first.")
            return

        print(f"\n{'='*100}")
        print("GENE ESSENTIALITY EVALUATION SUMMARY")
        print(f"{'='*100}")

        # Header
        print(f"{'Method':<35} {'Accuracy':>10} {'Precision':>10} {'Recall':>10} "
              f"{'F1-Score':>10} {'MCC':>10}")
        print("-" * 100)

        # Results
        for method_name, evaluation in self.essentiality_evaluations.items():
            print(f"{method_name:<35} {evaluation.accuracy:>10.4f} "
                  f"{evaluation.precision:>10.4f} {evaluation.recall:>10.4f} "
                  f"{evaluation.f1_score:>10.4f} {evaluation.mcc:>10.4f}")

        print("=" * 100)


def generate_sample_expression(model, condition: str = 'aerobic_glucose',
                               noise_level: float = 0.2) -> Dict[str, float]:
    """
    Generate sample gene expression data for benchmarking.

    Parameters
    ----------
    model : CBModel
        The metabolic model
    condition : str
        Growth condition: 'aerobic_glucose', 'anaerobic_glucose', 'aerobic_acetate'
    noise_level : float
        Standard deviation of random noise

    Returns
    -------
    dict
        Gene expression data
    """
    np.random.seed(42)
    base_expression = 5.0
    gene_exp = {}

    for gene_id in model.genes:
        expr = base_expression + np.random.normal(0, noise_level * base_expression)
        expr = max(0.1, expr)

        gene_lower = gene_id.lower()

        if condition == 'aerobic_glucose':
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh', 'pyk']):
                expr *= 3.0
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh', 'acn']):
                expr *= 2.5
            elif any(x in gene_lower for x in ['cyo', 'nuo', 'ndh']):
                expr *= 2.5
            elif any(x in gene_lower for x in ['pta', 'ack', 'adh', 'ldh']):
                expr *= 0.2

        elif condition == 'anaerobic_glucose':
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh', 'pyk']):
                expr *= 3.5
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh', 'acn']):
                expr *= 0.8
            elif any(x in gene_lower for x in ['cyo', 'nuo', 'ndh']):
                expr *= 0.1
            elif any(x in gene_lower for x in ['pta', 'ack', 'adh', 'ldh', 'pfl']):
                expr *= 4.0

        elif condition == 'aerobic_acetate':
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh']):
                expr *= 0.5
            elif any(x in gene_lower for x in ['pck', 'fbp', 'pps']):
                expr *= 3.0
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh']):
                expr *= 3.0
            elif any(x in gene_lower for x in ['aceb', 'acea']):
                expr *= 4.0

        gene_exp[gene_id] = expr

    return gene_exp


def main():
    """Main execution function demonstrating benchmark usage."""

    # Load model
    model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                             'tests', 'data')
    if os.path.exists('tests/data/e_coli_core.xml.gz'):
        model_path = 'tests/data/e_coli_core.xml.gz'
    else:
        model_path = os.path.join(model_dir, 'e_coli_core.xml.gz')

    print(f"Loading model: {model_path}")
    model = load_cbmodel(model_path)
    print(f"Model loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # Generate expression data
    print("\nGenerating sample gene expression data...")
    gene_exp = generate_sample_expression(model, 'aerobic_glucose')

    # ========================================================================
    # Example 1: Basic Benchmark
    # ========================================================================
    print("\n" + "="*70)
    print("EXAMPLE 1: Basic Method Comparison")
    print("="*70)

    benchmark1 = TranscriptomicsBenchmark(model, gene_exp)

    # Add methods to compare
    benchmark1.add_methods([
        FBAMethod(),
        pFBAMethod(),
        GIMMEMethod(cutoff=25, growth_frac=0.9),
        EFluxMethod(),
    ])

    # Run benchmark
    results1 = benchmark1.run_all(verbose=True)

    # Print summary
    benchmark1.print_summary()

    # Compare key fluxes
    key_reactions = ['R_GLCpts', 'R_PGI', 'R_PFK', 'R_PYK', 'R_CS',
                     'R_ACONTa', 'R_CYTBD', 'R_ACKr']
    benchmark1.compare_fluxes(key_reactions)

    # ========================================================================
    # Example 2: Parameter Sensitivity Analysis
    # ========================================================================
    print("\n\n" + "="*70)
    print("EXAMPLE 2: GIMME Parameter Sensitivity")
    print("="*70)

    benchmark2 = TranscriptomicsBenchmark(model, gene_exp)

    # Add GIMME with different parameters
    cutoffs = [10, 25, 50, 75]
    growth_fracs = [0.7, 0.9, 0.95]

    benchmark2.add_method(FBAMethod())  # Baseline

    for cutoff in cutoffs:
        benchmark2.add_method(GIMMEMethod(cutoff=cutoff, growth_frac=0.9))

    for gf in growth_fracs:
        benchmark2.add_method(GIMMEMethod(cutoff=25, growth_frac=gf))

    results2 = benchmark2.run_all(verbose=True)
    benchmark2.print_summary()

    # ========================================================================
    # Example 3: Multi-Condition Comparison
    # ========================================================================
    print("\n\n" + "="*70)
    print("EXAMPLE 3: Multi-Condition Analysis")
    print("="*70)

    conditions = ['aerobic_glucose', 'anaerobic_glucose', 'aerobic_acetate']

    for condition in conditions:
        print(f"\n--- Condition: {condition} ---")
        gene_exp_cond = generate_sample_expression(model, condition)

        benchmark3 = TranscriptomicsBenchmark(model, gene_exp_cond)
        benchmark3.add_methods([
            FBAMethod(),
            GIMMEMethod(cutoff=25, growth_frac=0.8),
            EFluxMethod(),
        ])

        benchmark3.run_all(verbose=False)
        benchmark3.print_summary()

    # ========================================================================
    # Example 4: Export Results
    # ========================================================================
    print("\n\n" + "="*70)
    print("EXAMPLE 4: Exporting Results")
    print("="*70)

    # Export summary
    benchmark1.export_results('benchmark_summary.csv', format='csv')

    # Export detailed fluxes
    benchmark1.export_detailed_fluxes('benchmark_fluxes.csv')

    # ========================================================================
    # Example 5: Top Varying Reactions
    # ========================================================================
    print("\n\n" + "="*70)
    print("EXAMPLE 5: Top Reactions by Variance")
    print("="*70)

    benchmark1.compare_fluxes([], top_n=20)

    print("\n" + "="*70)
    print("Benchmark Complete!")
    print("="*70)
    print("\nKey Takeaways:")
    print("1. Easy comparison of multiple methods with same data")
    print("2. Parameter sensitivity analysis for method optimization")
    print("3. Multi-condition benchmarking for robustness testing")
    print("4. Export capabilities for further analysis")
    print("5. Extensible framework for adding custom methods")


if __name__ == '__main__':
    main()
