#!/usr/bin/env python
"""
GIMME and E-Flux Example for E. coli Metabolic Models

This example demonstrates how to integrate transcriptomics data with
genome-scale metabolic models using GIMME and E-Flux methods.

Methods:
- GIMME: Gene Inactivity Moderated by Metabolism and Expression (Becker & Palsson, 2008)
- E-Flux: Expression-based Flux predictions (Colijn et al., 2009)

Author: ReFramed
Date: 2025
"""

import numpy as np
from reframed import load_cbmodel, FBA, GIMME, eFlux
from reframed.cobra.variability import FVA
import os

def generate_sample_expression_data(model, expression_level='high', noise_level=0.2):
    """
    Generate sample gene expression data for demonstration.

    In practice, you would load this from RNA-seq or microarray data.

    Parameters
    ----------
    model : CBModel
        The metabolic model
    expression_level : str
        'high', 'medium', or 'low' - overall expression level
    noise_level : float
        Standard deviation of random noise (default: 0.2)

    Returns
    -------
    dict
        Dictionary mapping gene IDs to expression values
    """
    np.random.seed(42)  # For reproducibility

    # Set base expression level
    if expression_level == 'high':
        base = 10.0
    elif expression_level == 'medium':
        base = 5.0
    else:  # low
        base = 2.0

    # Generate expression data with some realistic patterns
    gene_exp = {}
    for gene_id in model.genes:
        # Add random noise
        expression = base + np.random.normal(0, noise_level * base)

        # Ensure non-negative expression
        expression = max(0.1, expression)

        # Add some biological variation:
        # - Some genes highly expressed (e.g., glycolysis)
        # - Some genes lowly expressed (e.g., stress response)
        if 'pgi' in gene_id.lower() or 'pfk' in gene_id.lower() or 'gapdh' in gene_id.lower():
            expression *= 2.0  # Highly expressed glycolytic genes
        elif 'nuo' in gene_id.lower() or 'cyo' in gene_id.lower():
            expression *= 1.5  # Moderately expressed respiratory genes
        elif 'pta' in gene_id.lower() or 'ack' in gene_id.lower():
            expression *= 0.3  # Low expression of fermentation pathways

        gene_exp[gene_id] = expression

    return gene_exp


def simulate_condition_specific_expression(model, condition='aerobic_glucose'):
    """
    Generate condition-specific gene expression data.

    This simulates realistic expression patterns for different growth conditions.

    Parameters
    ----------
    model : CBModel
        The metabolic model
    condition : str
        Growth condition: 'aerobic_glucose', 'anaerobic_glucose', 'aerobic_acetate'

    Returns
    -------
    dict
        Gene expression dictionary
    """
    np.random.seed(42)
    base_expression = 5.0
    gene_exp = {}

    for gene_id in model.genes:
        expr = base_expression + np.random.normal(0, 1.0)
        expr = max(0.1, expr)

        gene_lower = gene_id.lower()

        if condition == 'aerobic_glucose':
            # High glycolysis and TCA cycle, high respiration
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh', 'pyk']):
                expr *= 3.0  # High glycolysis
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh', 'acn']):
                expr *= 2.5  # High TCA cycle
            elif any(x in gene_lower for x in ['cyo', 'nuo', 'ndh']):
                expr *= 2.5  # High respiratory chain
            elif any(x in gene_lower for x in ['pta', 'ack', 'adh', 'ldh']):
                expr *= 0.2  # Low fermentation

        elif condition == 'anaerobic_glucose':
            # High glycolysis, low TCA, no respiration, high fermentation
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh', 'pyk']):
                expr *= 3.5  # Very high glycolysis
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh', 'acn']):
                expr *= 0.8  # Reduced TCA cycle
            elif any(x in gene_lower for x in ['cyo', 'nuo', 'ndh']):
                expr *= 0.1  # Very low respiration
            elif any(x in gene_lower for x in ['pta', 'ack', 'adh', 'ldh', 'pfl']):
                expr *= 4.0  # High fermentation

        elif condition == 'aerobic_acetate':
            # Low glycolysis, high gluconeogenesis, high TCA and glyoxylate shunt
            if any(x in gene_lower for x in ['pgi', 'pfk', 'gapdh']):
                expr *= 0.5  # Low glycolysis
            elif any(x in gene_lower for x in ['pck', 'fbp', 'pps']):
                expr *= 3.0  # High gluconeogenesis
            elif any(x in gene_lower for x in ['mdh', 'fum', 'sdh']):
                expr *= 3.0  # High TCA cycle
            elif any(x in gene_lower for x in ['aceb', 'acea']):  # Glyoxylate shunt
                expr *= 4.0
            elif any(x in gene_lower for x in ['cyo', 'nuo']):
                expr *= 2.0  # Moderate respiration

        gene_exp[gene_id] = expr

    return gene_exp


def run_gimme_analysis(model, gene_exp, cutoff=25, growth_frac=0.9):
    """
    Run GIMME (Gene Inactivity Moderated by Metabolism and Expression) analysis.

    GIMME minimizes the use of reactions associated with lowly-expressed genes
    while maintaining a minimum growth rate.

    Parameters
    ----------
    model : CBModel
        Metabolic model
    gene_exp : dict
        Gene expression data {gene_id: expression_value}
    cutoff : float
        Percentile cutoff for expression threshold (default: 25)
    growth_frac : float
        Minimum growth as fraction of maximum (default: 0.9)

    Returns
    -------
    Solution
        GIMME solution object
    """
    print(f"\n{'='*70}")
    print("GIMME Analysis")
    print(f"{'='*70}")
    print(f"Number of genes with expression data: {len(gene_exp)}")
    print(f"Expression cutoff percentile: {cutoff}")
    print(f"Minimum growth fraction: {growth_frac}")

    # Run GIMME
    solution = GIMME(model, gene_exp, cutoff=cutoff, growth_frac=growth_frac)

    if solution.status.name == 'OPTIMAL':
        print(f"\nStatus: {solution.status.name}")
        print(f"Growth rate: {solution.fobj:.4f}")
        print(f"Biomass flux: {solution.values.get(model.biomass_reaction, 0):.4f}")

        # Show some key fluxes
        print("\nKey metabolic fluxes:")
        key_reactions = ['R_GLCpts', 'R_PGI', 'R_PFK', 'R_PYK', 'R_ACONTa',
                        'R_CS', 'R_CYTBD', 'R_ACKr', 'R_PTAr']
        for rxn_id in key_reactions:
            if rxn_id in solution.values:
                print(f"  {rxn_id:15s}: {solution.values[rxn_id]:8.4f}")
    else:
        print(f"GIMME failed with status: {solution.status.name}")

    return solution


def run_eflux_analysis(model, gene_exp, scale_rxn=None):
    """
    Run E-Flux (Expression-based Flux) analysis.

    E-Flux constrains reaction upper bounds proportionally to gene expression levels,
    then performs flux balance analysis.

    Parameters
    ----------
    model : CBModel
        Metabolic model
    gene_exp : dict
        Gene expression data {gene_id: expression_value}
    scale_rxn : str, optional
        Reaction ID to normalize fluxes to (e.g., glucose uptake)

    Returns
    -------
    Solution
        E-Flux solution object
    """
    print(f"\n{'='*70}")
    print("E-Flux Analysis")
    print(f"{'='*70}")
    print(f"Number of genes with expression data: {len(gene_exp)}")
    if scale_rxn:
        print(f"Scaling reaction: {scale_rxn}")

    # Run E-Flux
    solution = eFlux(model, gene_exp, scale_rxn=scale_rxn)

    if solution.status.name == 'OPTIMAL':
        print(f"\nStatus: {solution.status.name}")
        print(f"Growth rate: {solution.fobj:.4f}")
        print(f"Biomass flux: {solution.values.get(model.biomass_reaction, 0):.4f}")

        # Show some key fluxes
        print("\nKey metabolic fluxes:")
        key_reactions = ['R_GLCpts', 'R_PGI', 'R_PFK', 'R_PYK', 'R_ACONTa',
                        'R_CS', 'R_CYTBD', 'R_ACKr', 'R_PTAr']
        for rxn_id in key_reactions:
            if rxn_id in solution.values:
                print(f"  {rxn_id:15s}: {solution.values[rxn_id]:8.4f}")
    else:
        print(f"E-Flux failed with status: {solution.status.name}")

    return solution


def compare_methods(model, gene_exp, model_name="E. coli"):
    """
    Compare standard FBA, GIMME, and E-Flux results.

    Parameters
    ----------
    model : CBModel
        Metabolic model
    gene_exp : dict
        Gene expression data
    model_name : str
        Name of the model for display
    """
    print(f"\n{'#'*70}")
    print(f"Comparing FBA, GIMME, and E-Flux for {model_name}")
    print(f"{'#'*70}")

    # Standard FBA (baseline)
    print("\n[1] Standard FBA (no expression data)")
    print("-" * 70)
    fba_solution = FBA(model)
    print(f"Growth rate: {fba_solution.fobj:.4f}")

    # GIMME
    print("\n[2] GIMME (expression-based)")
    print("-" * 70)
    gimme_solution = run_gimme_analysis(model, gene_exp, cutoff=25, growth_frac=0.9)

    # E-Flux
    print("\n[3] E-Flux (expression-based)")
    print("-" * 70)
    eflux_solution = run_eflux_analysis(model, gene_exp)

    # Summary comparison
    print(f"\n{'='*70}")
    print("Summary Comparison")
    print(f"{'='*70}")
    print(f"{'Method':<20} {'Growth Rate':>15} {'Active Reactions':>20}")
    print("-" * 70)

    fba_active = sum(1 for v in fba_solution.values.values() if abs(v) > 1e-6)
    gimme_active = sum(1 for v in gimme_solution.values.values() if abs(v) > 1e-6)
    eflux_active = sum(1 for v in eflux_solution.values.values() if abs(v) > 1e-6)

    print(f"{'FBA':<20} {fba_solution.fobj:>15.4f} {fba_active:>20}")
    print(f"{'GIMME':<20} {gimme_solution.fobj:>15.4f} {gimme_active:>20}")
    print(f"{'E-Flux':<20} {eflux_solution.fobj:>15.4f} {eflux_active:>20}")

    return fba_solution, gimme_solution, eflux_solution


def analyze_flux_differences(fba_sol, gimme_sol, eflux_sol, model, top_n=10):
    """
    Analyze and display differences in flux predictions between methods.

    Parameters
    ----------
    fba_sol : Solution
        FBA solution
    gimme_sol : Solution
        GIMME solution
    eflux_sol : Solution
        E-Flux solution
    model : CBModel
        Metabolic model
    top_n : int
        Number of top differences to show
    """
    print(f"\n{'='*70}")
    print(f"Top {top_n} Flux Differences Between Methods")
    print(f"{'='*70}")

    # Calculate differences
    differences = []
    all_reactions = set(fba_sol.values.keys()) | set(gimme_sol.values.keys()) | set(eflux_sol.values.keys())

    for rxn_id in all_reactions:
        fba_flux = fba_sol.values.get(rxn_id, 0)
        gimme_flux = gimme_sol.values.get(rxn_id, 0)
        eflux_flux = eflux_sol.values.get(rxn_id, 0)

        # Calculate total variation
        variation = abs(fba_flux - gimme_flux) + abs(fba_flux - eflux_flux) + abs(gimme_flux - eflux_flux)

        if variation > 1e-3:  # Only consider significant differences
            differences.append((rxn_id, fba_flux, gimme_flux, eflux_flux, variation))

    # Sort by variation
    differences.sort(key=lambda x: x[4], reverse=True)

    # Display top differences
    print(f"\n{'Reaction ID':<15} {'FBA':>10} {'GIMME':>10} {'E-Flux':>10} {'Variation':>12}")
    print("-" * 70)
    for rxn_id, fba_f, gimme_f, eflux_f, var in differences[:top_n]:
        print(f"{rxn_id:<15} {fba_f:>10.4f} {gimme_f:>10.4f} {eflux_f:>10.4f} {var:>12.4f}")


def main():
    """Main execution function."""

    # Define model path
    model_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tests', 'data')

    # Check if running from examples/ or root directory
    if os.path.exists('tests/data/e_coli_core.xml.gz'):
        core_model_path = 'tests/data/e_coli_core.xml.gz'
    else:
        core_model_path = os.path.join(model_dir, 'e_coli_core.xml.gz')

    # ========================================================================
    # Example 1: E. coli Core Model
    # ========================================================================
    print("\n" + "="*70)
    print("Example 1: E. coli Core Model Analysis")
    print("="*70)

    print(f"\nLoading model: {core_model_path}")

    core_model = load_cbmodel(core_model_path)
    print(f"Model loaded: {len(core_model.reactions)} reactions, {len(core_model.metabolites)} metabolites, {len(core_model.genes)} genes")

    # Generate sample expression data for aerobic glucose condition
    print("\nGenerating sample gene expression data (aerobic glucose)...")
    gene_exp_aerobic = simulate_condition_specific_expression(core_model, 'aerobic_glucose')

    print(f"Expression data statistics:")
    exp_values = list(gene_exp_aerobic.values())
    print(f"  Mean: {np.mean(exp_values):.2f}")
    print(f"  Std:  {np.std(exp_values):.2f}")
    print(f"  Min:  {np.min(exp_values):.2f}")
    print(f"  Max:  {np.max(exp_values):.2f}")

    # Run comparison
    fba_sol, gimme_sol, eflux_sol = compare_methods(core_model, gene_exp_aerobic, "E. coli core")

    # Analyze differences
    analyze_flux_differences(fba_sol, gimme_sol, eflux_sol, core_model, top_n=15)

    # ========================================================================
    # Example 2: Condition-Specific Analysis (Anaerobic vs Aerobic)
    # ========================================================================
    print("\n\n" + "="*70)
    print("Example 2: Condition-Specific Analysis")
    print("="*70)

    # Anaerobic glucose
    print("\n--- Anaerobic Glucose Condition ---")
    gene_exp_anaerobic = simulate_condition_specific_expression(core_model, 'anaerobic_glucose')

    print("\nRunning GIMME for anaerobic condition...")
    gimme_anaerobic = run_gimme_analysis(core_model, gene_exp_anaerobic, cutoff=25, growth_frac=0.8)

    print("\nRunning E-Flux for anaerobic condition...")
    eflux_anaerobic = run_eflux_analysis(core_model, gene_exp_anaerobic)

    # Aerobic acetate
    print("\n--- Aerobic Acetate Condition ---")
    gene_exp_acetate = simulate_condition_specific_expression(core_model, 'aerobic_acetate')

    print("\nRunning GIMME for aerobic acetate condition...")
    gimme_acetate = run_gimme_analysis(core_model, gene_exp_acetate, cutoff=25, growth_frac=0.7)

    print("\nRunning E-Flux for aerobic acetate condition...")
    eflux_acetate = run_eflux_analysis(core_model, gene_exp_acetate)

    # ========================================================================
    # Example 3: Parameter Sensitivity Analysis
    # ========================================================================
    print("\n\n" + "="*70)
    print("Example 3: GIMME Parameter Sensitivity")
    print("="*70)

    print("\nTesting different cutoff percentiles...")
    cutoffs = [10, 25, 50, 75, 90]

    print(f"\n{'Cutoff %':<12} {'Growth Rate':>15} {'Active Rxns':>15}")
    print("-" * 50)

    for cutoff in cutoffs:
        sol = GIMME(core_model, gene_exp_aerobic, cutoff=cutoff, growth_frac=0.9)
        active = sum(1 for v in sol.values.values() if abs(v) > 1e-6)
        print(f"{cutoff:<12} {sol.fobj:>15.4f} {active:>15}")

    print("\nTesting different growth fraction requirements...")
    growth_fracs = [0.5, 0.7, 0.9, 0.95, 0.99]

    print(f"\n{'Growth Frac':<12} {'Growth Rate':>15} {'Active Rxns':>15}")
    print("-" * 50)

    for gf in growth_fracs:
        sol = GIMME(core_model, gene_exp_aerobic, cutoff=25, growth_frac=gf)
        active = sum(1 for v in sol.values.values() if abs(v) > 1e-6)
        print(f"{gf:<12.2f} {sol.fobj:>15.4f} {active:>15}")

    # ========================================================================
    # Example 4: Gene Expression Integration Workflow
    # ========================================================================
    print("\n\n" + "="*70)
    print("Example 4: Complete Workflow with Expression Data")
    print("="*70)

    print("\nStep 1: Load your expression data")
    print("  In practice, you would load from file:")
    print("    gene_exp = load_expression_from_file('expression_data.csv')")
    print("  Format: {'gene_id': expression_value, ...}")

    print("\nStep 2: Preprocess and quality check")
    print(f"  Total genes in model: {len(core_model.genes)}")
    print(f"  Genes with expression: {len(gene_exp_aerobic)}")
    print(f"  Coverage: {100*len(gene_exp_aerobic)/len(core_model.genes):.1f}%")

    print("\nStep 3: Run expression-integrated methods")
    print("  - GIMME: penalizes low-expression reactions")
    print("  - E-Flux: constrains bounds by expression")

    print("\nStep 4: Analyze and validate results")
    print("  - Compare to experimental flux measurements")
    print("  - Perform flux variability analysis")
    print("  - Test gene knockout predictions")

    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)

    print("\nKey Takeaways:")
    print("1. GIMME maintains minimum growth while minimizing low-expression fluxes")
    print("2. E-Flux directly constrains reaction bounds by expression levels")
    print("3. Both methods can capture condition-specific metabolism")
    print("4. Results should be validated with experimental data")

    print("\nFor more information:")
    print("- GIMME: Becker & Palsson (2008) PLoS Comput Biol")
    print("- E-Flux: Colijn et al. (2009) PLoS Comput Biol")
    print("- Documentation: https://reframed.readthedocs.io")


if __name__ == '__main__':
    main()
