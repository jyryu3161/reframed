#!/usr/bin/env python
"""
Template for Creating Custom Transcriptomics Integration Methods

이 파일은 새로운 전사체 통합 방법론을 구현하기 위한 템플릿입니다.
This file provides a template for implementing new transcriptomics integration methods.

Usage:
1. Copy this file and rename it (e.g., my_new_method.py)
2. Implement the CustomMethod class with your algorithm
3. Add your method to the benchmark
4. Run and compare results

Author: ReFramed
Date: 2025
"""

import time
import numpy as np
from typing import Dict, Optional, List, Tuple
from copy import deepcopy

from reframed import FBA
from reframed.cobra.transcriptomics import gene_to_reaction_expression
from reframed.solvers import solver_instance
from reframed.solvers.solution import Status

# Import benchmark framework
from benchmark_transcriptomics import TranscriptomicsMethod, BenchmarkResult


# ============================================================================
# TEMPLATE 1: Simple Expression-Based Method
# ============================================================================

class SimpleExpressionMethod(TranscriptomicsMethod):
    """
    Template for simple expression-based flux prediction methods.

    This template demonstrates how to:
    - Convert gene expression to reaction expression
    - Modify model constraints based on expression
    - Run FBA and return results

    Example: A method that sets lower bounds proportional to expression
    """

    def __init__(self, threshold: float = 1.0, scaling_factor: float = 0.1):
        """
        Initialize the method.

        Parameters
        ----------
        threshold : float
            Expression threshold below which reactions are constrained
        scaling_factor : float
            Scaling factor for expression-based constraints
        """
        super().__init__(f"SimpleExpression(t={threshold})")
        self.threshold = threshold
        self.scaling_factor = scaling_factor

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        """
        Execute the method.

        Step-by-step implementation:
        1. Convert gene expression to reaction expression
        2. Modify reaction bounds based on expression
        3. Run FBA
        4. Return results
        """
        try:
            start_time = time.time()

            # Step 1: Convert gene expression to reaction expression
            rxn_exp = gene_to_reaction_expression(
                model,
                gene_exp,
                and_func=min,  # For protein subunits (AND logic)
                or_func=max    # For isozymes (OR logic)
            )

            # Step 2: Create modified constraints
            # Example: Set lower bound based on expression level
            constraints = {}

            for rxn_id, expr_val in rxn_exp.items():
                if expr_val < self.threshold:
                    # Low expression: constrain the reaction
                    rxn = model.reactions[rxn_id]
                    new_bound = expr_val * self.scaling_factor

                    # Apply constraint based on reaction reversibility
                    if rxn.lb < 0:  # Reversible
                        constraints[rxn_id] = (-new_bound, new_bound)
                    else:  # Irreversible
                        constraints[rxn_id] = (0, new_bound)

            # Step 3: Run FBA with constraints
            solution = FBA(model, constraints=constraints)

            # Step 4: Calculate metrics and return
            exec_time = time.time() - start_time

            params = {
                'threshold': self.threshold,
                'scaling_factor': self.scaling_factor
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
                parameters={'threshold': self.threshold},
                error=str(e)
            )


# ============================================================================
# TEMPLATE 2: Optimization-Based Method
# ============================================================================

class OptimizationBasedMethod(TranscriptomicsMethod):
    """
    Template for optimization-based methods using custom objective functions.

    This template demonstrates how to:
    - Create a custom solver instance
    - Add custom variables and constraints
    - Define custom objective function
    - Solve and extract results

    Example: Minimize deviation from expression-predicted fluxes
    """

    def __init__(self, weight: float = 1.0, min_growth: float = 0.1):
        """
        Initialize the method.

        Parameters
        ----------
        weight : float
            Weight for expression deviation in objective
        min_growth : float
            Minimum required growth rate
        """
        super().__init__(f"OptimizationBased(w={weight})")
        self.weight = weight
        self.min_growth = min_growth

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        """
        Execute optimization-based method.

        Implementation steps:
        1. Convert expression to expected fluxes
        2. Create solver with custom objective
        3. Add constraints for minimum growth
        4. Solve and return results
        """
        try:
            start_time = time.time()

            # Step 1: Convert expression to expected flux levels
            rxn_exp = gene_to_reaction_expression(model, gene_exp)

            # Normalize expression to [0, 1]
            max_exp = max(rxn_exp.values()) if rxn_exp else 1.0
            expected_fluxes = {
                rxn_id: expr / max_exp
                for rxn_id, expr in rxn_exp.items()
            }

            # Step 2: Create solver
            solver = solver_instance(model)

            # Step 3: Set up optimization problem
            # Minimize: deviation from expected fluxes
            # Subject to: minimum growth constraint

            objective = {}

            for rxn_id in model.reactions:
                expected = expected_fluxes.get(rxn_id, 0.5)  # Default medium expression
                # Objective: minimize |flux - expected|
                # This is approximated as minimizing flux variance
                objective[rxn_id] = -expected * self.weight

            # Add minimum growth constraint
            if model.biomass_reaction:
                solver.add_constraint(
                    f'min_growth',
                    {model.biomass_reaction: 1},
                    '>', self.min_growth
                )

            # Step 4: Solve
            solution = solver.solve(
                objective,
                minimize=False,  # Maximize weighted flux
                get_values=True
            )

            exec_time = time.time() - start_time

            params = {
                'weight': self.weight,
                'min_growth': self.min_growth
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
                parameters={'weight': self.weight},
                error=str(e)
            )


# ============================================================================
# TEMPLATE 3: Two-Phase Method
# ============================================================================

class TwoPhaseMethod(TranscriptomicsMethod):
    """
    Template for two-phase optimization methods.

    This template demonstrates:
    - Phase 1: Optimize for growth
    - Phase 2: Optimize for expression consistency

    Example: Similar to GIMME but with custom phases
    """

    def __init__(self, expression_weight: float = 1.0, growth_frac: float = 0.9):
        """
        Initialize two-phase method.

        Parameters
        ----------
        expression_weight : float
            Weight for expression-based objective in phase 2
        growth_frac : float
            Minimum growth fraction (relative to max)
        """
        super().__init__(f"TwoPhase(w={expression_weight})")
        self.expression_weight = expression_weight
        self.growth_frac = growth_frac

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        """
        Execute two-phase optimization.

        Phase 1: Find maximum growth
        Phase 2: Optimize expression consistency with minimum growth constraint
        """
        try:
            start_time = time.time()

            # Phase 1: Maximum growth
            fba_solution = FBA(model)

            if fba_solution.status != Status.OPTIMAL:
                raise ValueError(f"Phase 1 failed: {fba_solution.status}")

            max_growth = fba_solution.fobj
            min_growth_required = max_growth * self.growth_frac

            # Phase 2: Optimize for expression consistency
            rxn_exp = gene_to_reaction_expression(model, gene_exp)

            # Calculate expression percentile threshold
            threshold = np.percentile(list(rxn_exp.values()), 25)

            # Create objective: penalize low-expression reactions
            objective = {}
            for rxn_id, expr_val in rxn_exp.items():
                if expr_val < threshold:
                    penalty = (threshold - expr_val) * self.expression_weight
                    objective[rxn_id] = penalty

            # Solve phase 2 with growth constraint
            constraints = {}
            if model.biomass_reaction:
                constraints[model.biomass_reaction] = (min_growth_required, float('inf'))

            solution = FBA(model, objective=objective, constraints=constraints, minimize=True)

            exec_time = time.time() - start_time

            params = {
                'expression_weight': self.expression_weight,
                'growth_frac': self.growth_frac,
                'max_growth': max_growth,
                'min_growth_required': min_growth_required
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
                parameters={'expression_weight': self.expression_weight},
                error=str(e)
            )


# ============================================================================
# TEMPLATE 4: Machine Learning-Based Method (Placeholder)
# ============================================================================

class MLBasedMethod(TranscriptomicsMethod):
    """
    Template for machine learning-based flux prediction.

    This template shows structure for ML-based methods that:
    - Use pre-trained models to predict fluxes
    - Integrate predictions as constraints
    - Solve modified FBA problem

    Note: This is a structural template. Implement ML model separately.
    """

    def __init__(self, model_path: Optional[str] = None, confidence_threshold: float = 0.7):
        """
        Initialize ML-based method.

        Parameters
        ----------
        model_path : str, optional
            Path to pre-trained ML model
        confidence_threshold : float
            Minimum confidence for using ML predictions
        """
        super().__init__(f"MLBased(conf={confidence_threshold})")
        self.model_path = model_path
        self.confidence_threshold = confidence_threshold
        self.ml_model = None  # Load your ML model here

    def run(self, model, gene_exp: Dict[str, float], **kwargs) -> BenchmarkResult:
        """
        Execute ML-based flux prediction.

        Steps:
        1. Preprocess expression data for ML model
        2. Predict flux ranges using ML model
        3. Use predictions as constraints in FBA
        4. Solve and return results
        """
        try:
            start_time = time.time()

            # Step 1: Preprocess expression data
            # TODO: Implement your preprocessing
            # features = self._preprocess_expression(gene_exp)

            # Step 2: Predict fluxes using ML model
            # TODO: Implement ML prediction
            # predictions = self.ml_model.predict(features)

            # Placeholder: Use expression as simple predictor
            rxn_exp = gene_to_reaction_expression(model, gene_exp)
            predictions = {rxn_id: expr for rxn_id, expr in rxn_exp.items()}

            # Step 3: Create constraints from predictions
            constraints = {}
            for rxn_id, predicted_flux in predictions.items():
                # Use predictions as soft constraints
                rxn = model.reactions[rxn_id]
                margin = predicted_flux * 0.2  # 20% margin

                if predicted_flux > 0:
                    constraints[rxn_id] = (
                        max(rxn.lb, predicted_flux - margin),
                        min(rxn.ub, predicted_flux + margin)
                    )

            # Step 4: Solve FBA with ML-based constraints
            solution = FBA(model, constraints=constraints)

            exec_time = time.time() - start_time

            params = {
                'model_path': self.model_path,
                'confidence_threshold': self.confidence_threshold
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
                parameters={'model_path': self.model_path},
                error=str(e)
            )

    def _preprocess_expression(self, gene_exp: Dict[str, float]) -> np.ndarray:
        """Preprocess expression data for ML model."""
        # TODO: Implement your preprocessing pipeline
        # Examples:
        # - Normalization
        # - Log transformation
        # - Feature engineering
        pass


# ============================================================================
# EXAMPLE: Using Custom Methods in Benchmark
# ============================================================================

def example_custom_method_usage():
    """
    Demonstrate how to use custom methods in the benchmark framework.
    """
    import os
    from reframed import load_cbmodel
    from benchmark_transcriptomics import (
        TranscriptomicsBenchmark,
        FBAMethod,
        GIMMEMethod,
        EFluxMethod,
        generate_sample_expression
    )

    # Load model
    if os.path.exists('tests/data/e_coli_core.xml.gz'):
        model_path = 'tests/data/e_coli_core.xml.gz'
    else:
        model_path = os.path.join('..', 'tests', 'data', 'e_coli_core.xml.gz')

    model = load_cbmodel(model_path)
    gene_exp = generate_sample_expression(model, 'aerobic_glucose')

    print("="*70)
    print("Testing Custom Methods")
    print("="*70)

    # Create benchmark
    benchmark = TranscriptomicsBenchmark(model, gene_exp)

    # Add standard methods
    benchmark.add_methods([
        FBAMethod(),
        GIMMEMethod(cutoff=25, growth_frac=0.9),
        EFluxMethod(),
    ])

    # Add custom methods
    benchmark.add_methods([
        SimpleExpressionMethod(threshold=2.0, scaling_factor=0.5),
        OptimizationBasedMethod(weight=1.0, min_growth=0.1),
        TwoPhaseMethod(expression_weight=1.0, growth_frac=0.9),
        # MLBasedMethod(confidence_threshold=0.7),  # Uncomment when ML model ready
    ])

    # Run benchmark
    results = benchmark.run_all(verbose=True)

    # Show results
    benchmark.print_summary()

    # Compare fluxes
    key_reactions = ['R_PGI', 'R_PFK', 'R_PYK', 'R_CS', 'R_CYTBD']
    benchmark.compare_fluxes(key_reactions)

    # Export results
    benchmark.export_results('custom_methods_benchmark.csv')

    print("\n" + "="*70)
    print("Custom Methods Test Complete!")
    print("="*70)


# ============================================================================
# CHECKLIST: Implementing Your Own Method
# ============================================================================

"""
□ Step 1: Choose a template that matches your method's approach
   - Simple constraint-based → Use Template 1
   - Custom optimization → Use Template 2
   - Multi-phase → Use Template 3
   - ML-based → Use Template 4

□ Step 2: Implement the __init__ method
   - Define all parameters your method needs
   - Set a descriptive name including key parameters
   - Store parameters as instance variables

□ Step 3: Implement the run method
   - Convert gene expression to reaction expression (if needed)
   - Implement your algorithm logic
   - Handle errors appropriately
   - Return BenchmarkResult with metrics

□ Step 4: Test your method
   - Run on small model first
   - Check for optimal solutions
   - Verify flux values make biological sense
   - Compare with baseline methods

□ Step 5: Add to benchmark
   - Import your method class
   - Create instance with parameters
   - Add to benchmark with add_method()
   - Run and analyze results

□ Step 6: Document your method
   - Add docstring explaining algorithm
   - Document all parameters
   - Include references to papers/methods
   - Add usage examples
"""


if __name__ == '__main__':
    example_custom_method_usage()
