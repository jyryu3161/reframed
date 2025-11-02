#!/usr/bin/env python
"""
Gene Essentiality Evaluation Module

This module provides tools to evaluate metabolic models using gene essentiality data
from knockout experiments (e.g., Keio collection for E. coli).

Features:
- Load essentiality data from CSV/JSON files
- Simulate single gene knockouts
- Compare predictions with experimental data
- Calculate performance metrics (accuracy, precision, recall, F1, MCC)
- Integrate with benchmark framework

Author: ReFramed
Date: 2025
"""

import os
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict

from reframed import FBA, pFBA, GIMME, eFlux
from reframed.cobra.knockout import gene_knockout


@dataclass
class GeneEssentiality:
    """Store gene essentiality information."""
    gene_id: str
    essentiality: str  # 'essential', 'non-essential', 'conditionally-essential'
    confidence: str = 'high'  # 'high', 'medium', 'low'
    condition: str = 'minimal_glucose'
    growth_rate_ratio: Optional[float] = None
    pathway: Optional[str] = None
    reference: Optional[str] = None


@dataclass
class EssentialityPrediction:
    """Store essentiality prediction result."""
    gene_id: str
    predicted_essential: bool
    growth_rate_wt: float
    growth_rate_ko: float
    growth_ratio: float
    status: str


@dataclass
class EssentialityEvaluation:
    """Store evaluation metrics."""
    method_name: str
    accuracy: float
    precision: float
    recall: float
    f1_score: float
    mcc: float  # Matthews Correlation Coefficient
    true_positives: int
    true_negatives: int
    false_positives: int
    false_negatives: int
    predictions: Dict[str, EssentialityPrediction] = field(default_factory=dict)
    threshold: float = 0.01


class EssentialityDataLoader:
    """Load and manage gene essentiality data."""

    @staticmethod
    def load_from_csv(filepath: str, condition: Optional[str] = None) -> Dict[str, GeneEssentiality]:
        """
        Load essentiality data from CSV file.

        Parameters
        ----------
        filepath : str
            Path to CSV file
        condition : str, optional
            Filter by growth condition

        Returns
        -------
        dict
            {gene_id: GeneEssentiality}
        """
        df = pd.read_csv(filepath, comment='#')

        essentiality_data = {}
        for _, row in df.iterrows():
            if condition and row['condition'] != condition:
                continue

            gene_ess = GeneEssentiality(
                gene_id=row['gene_id'],
                essentiality=row['essentiality'],
                confidence=row['confidence'],
                condition=row['condition'],
                reference=row.get('reference')
            )
            essentiality_data[row['gene_id']] = gene_ess

        return essentiality_data

    @staticmethod
    def load_from_json(filepath: str) -> Dict[str, GeneEssentiality]:
        """
        Load essentiality data from JSON file.

        Parameters
        ----------
        filepath : str
            Path to JSON file

        Returns
        -------
        dict
            {gene_id: GeneEssentiality}
        """
        with open(filepath, 'r') as f:
            data = json.load(f)

        essentiality_data = {}
        genes_data = data.get('genes', {})

        for gene_id, info in genes_data.items():
            gene_ess = GeneEssentiality(
                gene_id=gene_id,
                essentiality=info['essentiality'],
                confidence=info.get('confidence', 'high'),
                condition=data['metadata'].get('condition', 'unknown'),
                growth_rate_ratio=info.get('growth_rate_ratio'),
                pathway=info.get('pathway'),
                reference=data['metadata'].get('references', [None])[0]
            )
            essentiality_data[gene_id] = gene_ess

        return essentiality_data


class GeneEssentialityEvaluator:
    """
    Evaluate metabolic models using gene essentiality data.

    Example usage:
    >>> evaluator = GeneEssentialityEvaluator(model, essentiality_data)
    >>> evaluation = evaluator.evaluate_fba(threshold=0.01)
    >>> print(f"Accuracy: {evaluation.accuracy:.3f}")
    >>> evaluator.print_confusion_matrix(evaluation)
    """

    def __init__(self, model, essentiality_data: Dict[str, GeneEssentiality],
                 growth_threshold: float = 0.01):
        """
        Initialize evaluator.

        Parameters
        ----------
        model : CBModel
            Metabolic model
        essentiality_data : dict
            Gene essentiality data {gene_id: GeneEssentiality}
        growth_threshold : float
            Growth rate threshold for essentiality (default: 0.01)
            Genes with knockout growth < threshold are predicted essential
        """
        self.model = model
        self.essentiality_data = essentiality_data
        self.growth_threshold = growth_threshold

        # Filter for genes in both model and data
        self.evaluated_genes = set(essentiality_data.keys()) & set(model.genes.keys())

        if not self.evaluated_genes:
            raise ValueError("No genes in common between model and essentiality data")

    def predict_essentiality(self, method: str = 'fba', gene_exp: Optional[Dict] = None,
                            **method_kwargs) -> Dict[str, EssentialityPrediction]:
        """
        Predict gene essentiality using specified method.

        Parameters
        ----------
        method : str
            Method to use: 'fba', 'pfba', 'gimme', 'eflux'
        gene_exp : dict, optional
            Gene expression data (required for gimme/eflux)
        **method_kwargs
            Additional method-specific parameters

        Returns
        -------
        dict
            {gene_id: EssentialityPrediction}
        """
        # Get wild-type growth
        if method == 'fba':
            wt_solution = FBA(self.model)
        elif method == 'pfba':
            wt_solution = pFBA(self.model)
        elif method == 'gimme':
            if gene_exp is None:
                raise ValueError("gene_exp required for GIMME")
            wt_solution = GIMME(self.model, gene_exp, **method_kwargs)
        elif method == 'eflux':
            if gene_exp is None:
                raise ValueError("gene_exp required for E-Flux")
            wt_solution = eFlux(self.model, gene_exp, **method_kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")

        wt_growth = wt_solution.fobj

        # Predict essentiality for each gene
        predictions = {}

        for gene_id in self.evaluated_genes:
            # Simulate gene knockout
            try:
                if method == 'fba':
                    ko_solution = gene_knockout(self.model, [gene_id], method='FBA')
                elif method == 'pfba':
                    ko_solution = gene_knockout(self.model, [gene_id], method='pFBA')
                elif method == 'gimme':
                    # For GIMME/E-Flux, we need to do knockout then run the method
                    # Create a temporary knockout by setting gene expression to 0
                    ko_gene_exp = gene_exp.copy()
                    ko_gene_exp[gene_id] = 0.0
                    ko_solution = GIMME(self.model, ko_gene_exp, **method_kwargs)
                elif method == 'eflux':
                    ko_gene_exp = gene_exp.copy()
                    ko_gene_exp[gene_id] = 0.0
                    ko_solution = eFlux(self.model, ko_gene_exp, **method_kwargs)

                ko_growth = ko_solution.fobj
                growth_ratio = ko_growth / wt_growth if wt_growth > 1e-6 else 0

                predicted_essential = growth_ratio < self.growth_threshold

                predictions[gene_id] = EssentialityPrediction(
                    gene_id=gene_id,
                    predicted_essential=predicted_essential,
                    growth_rate_wt=wt_growth,
                    growth_rate_ko=ko_growth,
                    growth_ratio=growth_ratio,
                    status=ko_solution.status.name
                )

            except Exception as e:
                # If knockout fails, assume essential
                predictions[gene_id] = EssentialityPrediction(
                    gene_id=gene_id,
                    predicted_essential=True,
                    growth_rate_wt=wt_growth,
                    growth_rate_ko=0.0,
                    growth_ratio=0.0,
                    status=f"ERROR: {str(e)}"
                )

        return predictions

    def calculate_metrics(self, predictions: Dict[str, EssentialityPrediction],
                         method_name: str = "Unknown") -> EssentialityEvaluation:
        """
        Calculate evaluation metrics.

        Parameters
        ----------
        predictions : dict
            Prediction results {gene_id: EssentialityPrediction}
        method_name : str
            Name of the method

        Returns
        -------
        EssentialityEvaluation
            Evaluation metrics
        """
        tp = fp = tn = fn = 0

        for gene_id, pred in predictions.items():
            if gene_id not in self.essentiality_data:
                continue

            experimental = self.essentiality_data[gene_id]
            is_essential = experimental.essentiality == 'essential'

            if pred.predicted_essential and is_essential:
                tp += 1
            elif pred.predicted_essential and not is_essential:
                fp += 1
            elif not pred.predicted_essential and not is_essential:
                tn += 1
            elif not pred.predicted_essential and is_essential:
                fn += 1

        # Calculate metrics
        total = tp + tn + fp + fn
        accuracy = (tp + tn) / total if total > 0 else 0

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        # Matthews Correlation Coefficient
        numerator = (tp * tn) - (fp * fn)
        denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = numerator / denominator if denominator > 0 else 0

        return EssentialityEvaluation(
            method_name=method_name,
            accuracy=accuracy,
            precision=precision,
            recall=recall,
            f1_score=f1_score,
            mcc=mcc,
            true_positives=tp,
            true_negatives=tn,
            false_positives=fp,
            false_negatives=fn,
            predictions=predictions,
            threshold=self.growth_threshold
        )

    def evaluate_fba(self) -> EssentialityEvaluation:
        """Evaluate using FBA."""
        predictions = self.predict_essentiality(method='fba')
        return self.calculate_metrics(predictions, "FBA")

    def evaluate_pfba(self) -> EssentialityEvaluation:
        """Evaluate using pFBA."""
        predictions = self.predict_essentiality(method='pfba')
        return self.calculate_metrics(predictions, "pFBA")

    def evaluate_gimme(self, gene_exp: Dict, **kwargs) -> EssentialityEvaluation:
        """Evaluate using GIMME."""
        predictions = self.predict_essentiality(method='gimme', gene_exp=gene_exp, **kwargs)
        return self.calculate_metrics(predictions, f"GIMME({kwargs})")

    def evaluate_eflux(self, gene_exp: Dict, **kwargs) -> EssentialityEvaluation:
        """Evaluate using E-Flux."""
        predictions = self.predict_essentiality(method='eflux', gene_exp=gene_exp, **kwargs)
        return self.calculate_metrics(predictions, f"E-Flux({kwargs})")

    @staticmethod
    def print_confusion_matrix(evaluation: EssentialityEvaluation):
        """Print confusion matrix."""
        print(f"\n{'='*60}")
        print(f"Confusion Matrix: {evaluation.method_name}")
        print(f"{'='*60}")
        print(f"                  Predicted Essential  Predicted Non-Essential")
        print(f"Actual Essential        {evaluation.true_positives:5d}                 {evaluation.false_negatives:5d}")
        print(f"Actual Non-Essential    {evaluation.false_positives:5d}                 {evaluation.true_negatives:5d}")
        print(f"{'='*60}")

    @staticmethod
    def print_metrics(evaluation: EssentialityEvaluation):
        """Print evaluation metrics."""
        print(f"\n{'='*60}")
        print(f"Evaluation Metrics: {evaluation.method_name}")
        print(f"{'='*60}")
        print(f"Accuracy:   {evaluation.accuracy:.4f}  ({evaluation.accuracy*100:.2f}%)")
        print(f"Precision:  {evaluation.precision:.4f}  (TP / (TP + FP))")
        print(f"Recall:     {evaluation.recall:.4f}  (TP / (TP + FN))")
        print(f"F1-Score:   {evaluation.f1_score:.4f}")
        print(f"MCC:        {evaluation.mcc:.4f}  (Matthews Correlation)")
        print(f"{'='*60}")
        print(f"True Positives:   {evaluation.true_positives:4d}  (Correctly predicted essential)")
        print(f"True Negatives:   {evaluation.true_negatives:4d}  (Correctly predicted non-essential)")
        print(f"False Positives:  {evaluation.false_positives:4d}  (Incorrectly predicted essential)")
        print(f"False Negatives:  {evaluation.false_negatives:4d}  (Incorrectly predicted non-essential)")
        print(f"{'='*60}")

    def compare_methods(self, methods: List[str], gene_exp: Optional[Dict] = None,
                       method_params: Optional[Dict] = None) -> pd.DataFrame:
        """
        Compare multiple methods.

        Parameters
        ----------
        methods : list
            List of methods: ['fba', 'pfba', 'gimme', 'eflux']
        gene_exp : dict, optional
            Gene expression data (for gimme/eflux)
        method_params : dict, optional
            Method-specific parameters {method: kwargs}

        Returns
        -------
        DataFrame
            Comparison table
        """
        if method_params is None:
            method_params = {}

        results = []

        for method in methods:
            print(f"\nEvaluating {method.upper()}...")

            if method == 'fba':
                eval_result = self.evaluate_fba()
            elif method == 'pfba':
                eval_result = self.evaluate_pfba()
            elif method == 'gimme':
                params = method_params.get('gimme', {'cutoff': 25, 'growth_frac': 0.9})
                eval_result = self.evaluate_gimme(gene_exp, **params)
            elif method == 'eflux':
                params = method_params.get('eflux', {})
                eval_result = self.evaluate_eflux(gene_exp, **params)
            else:
                continue

            results.append({
                'Method': eval_result.method_name,
                'Accuracy': eval_result.accuracy,
                'Precision': eval_result.precision,
                'Recall': eval_result.recall,
                'F1-Score': eval_result.f1_score,
                'MCC': eval_result.mcc,
                'TP': eval_result.true_positives,
                'TN': eval_result.true_negatives,
                'FP': eval_result.false_positives,
                'FN': eval_result.false_negatives,
            })

        return pd.DataFrame(results)

    def analyze_misclassifications(self, evaluation: EssentialityEvaluation,
                                   max_display: int = 10):
        """
        Analyze and display misclassified genes.

        Parameters
        ----------
        evaluation : EssentialityEvaluation
            Evaluation results
        max_display : int
            Maximum number of genes to display per category
        """
        false_positives = []
        false_negatives = []

        for gene_id, pred in evaluation.predictions.items():
            if gene_id not in self.essentiality_data:
                continue

            experimental = self.essentiality_data[gene_id]
            is_essential = experimental.essentiality == 'essential'

            if pred.predicted_essential and not is_essential:
                false_positives.append((gene_id, pred, experimental))
            elif not pred.predicted_essential and is_essential:
                false_negatives.append((gene_id, pred, experimental))

        # Print false positives
        print(f"\n{'='*80}")
        print(f"False Positives: Predicted Essential but Actually Non-Essential")
        print(f"{'='*80}")
        print(f"{'Gene ID':<15} {'Growth Ratio':>12} {'Pathway':<20} {'Confidence':<12}")
        print("-" * 80)

        for gene_id, pred, exp in false_positives[:max_display]:
            pathway = exp.pathway or 'unknown'
            print(f"{gene_id:<15} {pred.growth_ratio:>12.4f} {pathway:<20} {exp.confidence:<12}")

        if len(false_positives) > max_display:
            print(f"... and {len(false_positives) - max_display} more")

        # Print false negatives
        print(f"\n{'='*80}")
        print(f"False Negatives: Predicted Non-Essential but Actually Essential")
        print(f"{'='*80}")
        print(f"{'Gene ID':<15} {'Growth Ratio':>12} {'Pathway':<20} {'Confidence':<12}")
        print("-" * 80)

        for gene_id, pred, exp in false_negatives[:max_display]:
            pathway = exp.pathway or 'unknown'
            print(f"{gene_id:<15} {pred.growth_ratio:>12.4f} {pathway:<20} {exp.confidence:<12}")

        if len(false_negatives) > max_display:
            print(f"... and {len(false_negatives) - max_display} more")

    def export_results(self, evaluation: EssentialityEvaluation, filename: str):
        """
        Export detailed results to CSV.

        Parameters
        ----------
        evaluation : EssentialityEvaluation
            Evaluation results
        filename : str
            Output filename
        """
        data = []

        for gene_id, pred in evaluation.predictions.items():
            if gene_id not in self.essentiality_data:
                continue

            experimental = self.essentiality_data[gene_id]

            data.append({
                'gene_id': gene_id,
                'experimental_essentiality': experimental.essentiality,
                'predicted_essential': pred.predicted_essential,
                'growth_rate_wt': pred.growth_rate_wt,
                'growth_rate_ko': pred.growth_rate_ko,
                'growth_ratio': pred.growth_ratio,
                'confidence': experimental.confidence,
                'pathway': experimental.pathway or '',
                'correct_prediction': (pred.predicted_essential == (experimental.essentiality == 'essential'))
            })

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)
        print(f"\nResults exported to: {filename}")


def load_essentiality_data(filepath: str, format: str = 'auto') -> Dict[str, GeneEssentiality]:
    """
    Load gene essentiality data from file.

    Parameters
    ----------
    filepath : str
        Path to data file
    format : str
        File format: 'csv', 'json', or 'auto' (detect from extension)

    Returns
    -------
    dict
        {gene_id: GeneEssentiality}
    """
    if format == 'auto':
        if filepath.endswith('.csv'):
            format = 'csv'
        elif filepath.endswith('.json'):
            format = 'json'
        else:
            raise ValueError(f"Cannot auto-detect format from: {filepath}")

    loader = EssentialityDataLoader()

    if format == 'csv':
        return loader.load_from_csv(filepath)
    elif format == 'json':
        return loader.load_from_json(filepath)
    else:
        raise ValueError(f"Unknown format: {format}")
