#!/usr/bin/env python
"""
Gene Essentiality Benchmark Example

Demonstrates how to evaluate metabolic models using gene essentiality data.

This example:
1. Loads essentiality data from experimental knockouts
2. Simulates gene knockouts using different methods
3. Compares predictions with experimental data
4. Calculates performance metrics

Author: ReFramed
Date: 2025
"""

import os
import sys
import numpy as np
from reframed import load_cbmodel

# Import essentiality evaluation module
from gene_essentiality_evaluation import (
    load_essentiality_data,
    GeneEssentialityEvaluator
)

# Import benchmark framework for expression data generation
from benchmark_transcriptomics import generate_sample_expression


def main():
    """Main execution function."""

    print("="*80)
    print("Gene Essentiality Evaluation Example")
    print("="*80)

    # ========================================================================
    # Step 1: Load Model
    # ========================================================================
    print("\nStep 1: Loading metabolic model...")

    if os.path.exists('tests/data/e_coli_core.xml.gz'):
        model_path = 'tests/data/e_coli_core.xml.gz'
    else:
        model_path = os.path.join('..', 'tests', 'data', 'e_coli_core.xml.gz')

    model = load_cbmodel(model_path)
    print(f"✓ Model loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # ========================================================================
    # Step 2: Load Essentiality Data
    # ========================================================================
    print("\nStep 2: Loading gene essentiality data...")

    # Try both CSV and JSON formats
    if os.path.exists('data/ecoli_gene_essentiality.csv'):
        ess_data_path = 'data/ecoli_gene_essentiality.csv'
    elif os.path.exists('examples/data/ecoli_gene_essentiality.csv'):
        ess_data_path = 'examples/data/ecoli_gene_essentiality.csv'
    else:
        print("✗ Error: Could not find essentiality data file")
        print("  Expected: data/ecoli_gene_essentiality.csv")
        return

    essentiality_data = load_essentiality_data(ess_data_path, format='csv')

    print(f"✓ Loaded essentiality data for {len(essentiality_data)} genes")

    # Count essential vs non-essential
    essential_count = sum(1 for e in essentiality_data.values() if e.essentiality == 'essential')
    non_essential_count = len(essentiality_data) - essential_count

    print(f"  - Essential genes: {essential_count}")
    print(f"  - Non-essential genes: {non_essential_count}")

    # ========================================================================
    # Step 3: Initialize Evaluator
    # ========================================================================
    print("\nStep 3: Initializing evaluator...")

    evaluator = GeneEssentialityEvaluator(
        model,
        essentiality_data,
        growth_threshold=0.01  # Genes with <1% growth are essential
    )

    print(f"✓ Evaluator ready")
    print(f"  - Genes to evaluate: {len(evaluator.evaluated_genes)}")
    print(f"  - Growth threshold: {evaluator.growth_threshold}")

    # ========================================================================
    # Example 1: Basic FBA Evaluation
    # ========================================================================
    print("\n" + "="*80)
    print("Example 1: FBA Evaluation")
    print("="*80)

    fba_eval = evaluator.evaluate_fba()

    evaluator.print_confusion_matrix(fba_eval)
    evaluator.print_metrics(fba_eval)

    # ========================================================================
    # Example 2: Compare Multiple Methods
    # ========================================================================
    print("\n\n" + "="*80)
    print("Example 2: Comparing Multiple Methods")
    print("="*80)

    # Generate sample expression data
    print("\nGenerating sample gene expression data...")
    gene_exp = generate_sample_expression(model, 'aerobic_glucose')
    print(f"✓ Expression data for {len(gene_exp)} genes")

    # Compare methods
    methods = ['fba', 'pfba']

    # For GIMME/E-Flux, we can add them if expression data is provided
    # Uncomment to include (warning: takes longer)
    # methods.extend(['gimme', 'eflux'])

    method_params = {
        'gimme': {'cutoff': 25, 'growth_frac': 0.9},
        'eflux': {}
    }

    comparison_df = evaluator.compare_methods(
        methods,
        gene_exp=gene_exp,
        method_params=method_params
    )

    print("\n" + "="*80)
    print("Method Comparison Summary")
    print("="*80)
    print(comparison_df.to_string(index=False))

    # ========================================================================
    # Example 3: Analyze Misclassifications
    # ========================================================================
    print("\n\n" + "="*80)
    print("Example 3: Analyzing Misclassifications")
    print("="*80)

    evaluator.analyze_misclassifications(fba_eval, max_display=10)

    # ========================================================================
    # Example 4: Threshold Sensitivity Analysis
    # ========================================================================
    print("\n\n" + "="*80)
    print("Example 4: Growth Threshold Sensitivity Analysis")
    print("="*80)

    thresholds = [0.001, 0.005, 0.01, 0.05, 0.1]
    threshold_results = []

    print("\nTesting different growth thresholds...")
    print(f"{'Threshold':>10} {'Accuracy':>10} {'Precision':>10} {'Recall':>10} "
          f"{'F1-Score':>10} {'MCC':>10}")
    print("-" * 70)

    for threshold in thresholds:
        evaluator.growth_threshold = threshold
        predictions = evaluator.predict_essentiality(method='fba')
        eval_result = evaluator.calculate_metrics(predictions, f"FBA(t={threshold})")

        print(f"{threshold:>10.3f} {eval_result.accuracy:>10.4f} "
              f"{eval_result.precision:>10.4f} {eval_result.recall:>10.4f} "
              f"{eval_result.f1_score:>10.4f} {eval_result.mcc:>10.4f}")

        threshold_results.append({
            'threshold': threshold,
            'accuracy': eval_result.accuracy,
            'f1_score': eval_result.f1_score,
            'mcc': eval_result.mcc
        })

    # Find best threshold
    best_threshold = max(threshold_results, key=lambda x: x['f1_score'])
    print(f"\n✓ Best threshold (by F1-score): {best_threshold['threshold']}")
    print(f"  Accuracy: {best_threshold['accuracy']:.4f}")
    print(f"  F1-Score: {best_threshold['f1_score']:.4f}")
    print(f"  MCC: {best_threshold['mcc']:.4f}")

    # Reset to default threshold
    evaluator.growth_threshold = 0.01

    # ========================================================================
    # Example 5: Export Results
    # ========================================================================
    print("\n\n" + "="*80)
    print("Example 5: Exporting Results")
    print("="*80)

    # Export detailed results
    evaluator.export_results(fba_eval, 'essentiality_predictions.csv')

    # Export comparison
    comparison_df.to_csv('essentiality_comparison.csv', index=False)
    print("Comparison table exported to: essentiality_comparison.csv")

    # ========================================================================
    # Example 6: Pathway-Specific Analysis
    # ========================================================================
    print("\n\n" + "="*80)
    print("Example 6: Pathway-Specific Analysis")
    print("="*80)

    # Load JSON data for pathway information
    if os.path.exists('data/ecoli_gene_essentiality.json'):
        json_path = 'data/ecoli_gene_essentiality.json'
    elif os.path.exists('examples/data/ecoli_gene_essentiality.json'):
        json_path = 'examples/data/ecoli_gene_essentiality.json'
    else:
        json_path = None

    if json_path:
        essentiality_data_json = load_essentiality_data(json_path, format='json')

        # Group by pathway
        from collections import defaultdict
        pathway_genes = defaultdict(list)

        for gene_id, gene_data in essentiality_data_json.items():
            if gene_data.pathway:
                pathway_genes[gene_data.pathway].append(gene_id)

        # Analyze each pathway
        print("\nEssentiality by Pathway:")
        print(f"{'Pathway':<30} {'Total':>8} {'Essential':>12} {'% Essential':>12}")
        print("-" * 70)

        for pathway, genes in sorted(pathway_genes.items()):
            essential = sum(1 for g in genes
                          if essentiality_data_json[g].essentiality == 'essential')
            pct = 100 * essential / len(genes) if genes else 0
            print(f"{pathway:<30} {len(genes):>8} {essential:>12} {pct:>11.1f}%")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n\n" + "="*80)
    print("✓ Essentiality Evaluation Complete!")
    print("="*80)

    print("\nKey Findings:")
    print(f"1. FBA Accuracy: {fba_eval.accuracy*100:.2f}%")
    print(f"2. Best performing method: {comparison_df.loc[comparison_df['Accuracy'].idxmax(), 'Method']}")
    print(f"3. Optimal threshold: {best_threshold['threshold']}")

    print("\nOutput Files:")
    print("- essentiality_predictions.csv: Detailed predictions for each gene")
    print("- essentiality_comparison.csv: Method comparison table")

    print("\nNext Steps:")
    print("1. Compare with different growth conditions")
    print("2. Integrate with transcriptomics data (GIMME/E-Flux)")
    print("3. Analyze pathway-specific essentiality patterns")
    print("4. Validate model predictions with experimental data")


if __name__ == '__main__':
    main()
