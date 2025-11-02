#!/usr/bin/env python
"""
Quick Benchmark Example

빠른 벤치마크 실행 예제 - 여러 방법론을 동시에 비교합니다.
Quick benchmark execution - Compare multiple methods simultaneously.

This script demonstrates the simplest way to benchmark transcriptomics
integration methods using the ReFramed framework.

Author: ReFramed
Date: 2025
"""

import os
from reframed import load_cbmodel
from benchmark_transcriptomics import (
    TranscriptomicsBenchmark,
    FBAMethod,
    pFBAMethod,
    GIMMEMethod,
    EFluxMethod,
    generate_sample_expression
)


def main():
    """Quick benchmark demonstration."""

    print("="*70)
    print("Quick Benchmark: Transcriptomics Integration Methods")
    print("="*70)

    # 1. Load model
    if os.path.exists('tests/data/e_coli_core.xml.gz'):
        model_path = 'tests/data/e_coli_core.xml.gz'
    else:
        model_path = os.path.join('..', 'tests', 'data', 'e_coli_core.xml.gz')

    print(f"\nLoading model: {model_path}")
    model = load_cbmodel(model_path)
    print(f"✓ Loaded: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, "
          f"{len(model.genes)} genes")

    # 2. Generate expression data
    print("\nGenerating sample gene expression data...")
    gene_exp = generate_sample_expression(model, 'aerobic_glucose')
    print(f"✓ Expression data for {len(gene_exp)} genes")

    # 3. Create benchmark
    benchmark = TranscriptomicsBenchmark(model, gene_exp)

    # 4. Add methods to compare
    print("\nAdding methods to benchmark:")
    methods = [
        ("FBA (baseline)", FBAMethod()),
        ("pFBA (parsimonious)", pFBAMethod()),
        ("GIMME (cutoff=25)", GIMMEMethod(cutoff=25, growth_frac=0.9)),
        ("GIMME (cutoff=50)", GIMMEMethod(cutoff=50, growth_frac=0.9)),
        ("E-Flux", EFluxMethod()),
    ]

    for desc, method in methods:
        print(f"  • {desc}")
        benchmark.add_method(method)

    # 5. Run benchmark
    print("\n" + "="*70)
    print("Running Benchmark...")
    print("="*70)
    results = benchmark.run_all(verbose=True)

    # 6. Show summary
    benchmark.print_summary()

    # 7. Compare key fluxes
    print("\n")
    key_reactions = ['R_GLCpts', 'R_PGI', 'R_PFK', 'R_PYK', 'R_CS',
                     'R_ACONTa', 'R_CYTBD', 'R_ACKr', 'R_PTAr']
    benchmark.compare_fluxes(key_reactions)

    # 8. Show most variable reactions
    print("\n")
    benchmark.compare_fluxes([], top_n=15)

    # 9. Export results
    print("\n" + "="*70)
    print("Exporting Results...")
    print("="*70)
    benchmark.export_results('benchmark_results.csv', format='csv')
    benchmark.export_detailed_fluxes('benchmark_detailed_fluxes.csv')

    print("\n" + "="*70)
    print("✓ Benchmark Complete!")
    print("="*70)
    print("\nNext Steps:")
    print("1. Review benchmark_results.csv for summary statistics")
    print("2. Analyze benchmark_detailed_fluxes.csv for flux comparisons")
    print("3. Add your own methods using custom_method_template.py")
    print("4. Run benchmark_transcriptomics.py for advanced examples")


if __name__ == '__main__':
    main()
