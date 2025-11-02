#!/usr/bin/env python
"""
Simple GIMME and E-Flux Example

Quick start guide for integrating gene expression with metabolic models.
"""

from reframed import load_cbmodel, FBA, GIMME, eFlux
import os

# Load E. coli core model
# Get the correct path whether running from examples/ or root directory
if os.path.exists('tests/data/e_coli_core.xml.gz'):
    model_path = 'tests/data/e_coli_core.xml.gz'
else:
    model_path = os.path.join('..', 'tests', 'data', 'e_coli_core.xml.gz')

model = load_cbmodel(model_path)

print(f"Loaded model with {len(model.genes)} genes")

# Example gene expression data
# In practice, load this from RNA-seq or microarray data
# Note: Gene IDs in the model have 'G_' prefix
gene_expression = {
    'G_b0008': 10.5,  # talB - transaldolase
    'G_b0114': 8.2,   # aceE - pyruvate dehydrogenase
    'G_b0115': 8.5,   # aceF - pyruvate dehydrogenase
    'G_b0116': 7.9,   # lpd - dihydrolipoamide dehydrogenase
    'G_b0118': 12.3,  # acnB - aconitate hydratase
    'G_b0356': 6.5,   # frmA - S-formylglutathione hydrolase
    'G_b0451': 15.2,  # cyoA - cytochrome o ubiquinol oxidase
    'G_b0720': 9.8,   # gltA - citrate synthase
    'G_b0721': 3.2,   # sdhC - succinate dehydrogenase
    'G_b0722': 3.5,   # sdhD - succinate dehydrogenase
    'G_b0723': 3.8,   # sdhA - succinate dehydrogenase
    'G_b0724': 3.1,   # sdhB - succinate dehydrogenase
    'G_b0726': 11.5,  # sucA - 2-oxoglutarate dehydrogenase
    'G_b0727': 11.2,  # sucB - dihydrolipoamide succinyltransferase
    'G_b0728': 10.8,  # sucC - succinyl-CoA synthetase
    'G_b0729': 10.5,  # sucD - succinyl-CoA synthetase
    'G_b0755': 4.2,   # galP - galactose permease
    'G_b0767': 2.1,   # maeB - malic enzyme
    'G_b0810': 1.8,   # tnaA - tryptophanase
    'G_b0875': 14.5,  # pykF - pyruvate kinase
    'G_b1136': 16.2,  # icd - isocitrate dehydrogenase
    'G_b1276': 13.8,  # acnA - aconitate hydratase
    'G_b1478': 5.5,   # adhE - alcohol dehydrogenase
    'G_b1602': 12.1,  # lpdA - dihydrolipoamide dehydrogenase
    'G_b1603': 11.8,  # frdA - fumarate reductase
    'G_b1604': 11.5,  # frdB - fumarate reductase
    'G_b1605': 11.2,  # frdC - fumarate reductase
    'G_b1676': 13.5,  # pykA - pyruvate kinase
    'G_b1779': 15.8,  # gapA - glyceraldehyde 3-phosphate dehydrogenase
    'G_b1817': 2.8,   # pflB - pyruvate formate-lyase
    'G_b1818': 2.5,   # focA - formate transporter
    'G_b1819': 14.2,  # pfkA - phosphofructokinase
    'G_b1852': 3.5,   # poxB - pyruvate oxidase
    'G_b2097': 7.8,   # mdh - malate dehydrogenase
    'G_b2114': 8.5,   # pfkB - phosphofructokinase
    'G_b2133': 16.5,  # gapC - glyceraldehyde 3-phosphate dehydrogenase
    'G_b2296': 5.2,   # ackA - acetate kinase
    'G_b2297': 5.0,   # pta - phosphate acetyltransferase
    'G_b2415': 2.2,   # cydA - cytochrome d ubiquinol oxidase
    'G_b2416': 2.1,   # cydB - cytochrome d ubiquinol oxidase
    'G_b2458': 9.5,   # fbaA - fructose-bisphosphate aldolase
    'G_b2465': 10.2,  # ppc - phosphoenolpyruvate carboxylase
    'G_b2779': 7.8,   # eno - enolase
    'G_b2925': 15.5,  # pgk - phosphoglycerate kinase
    'G_b2926': 14.8,  # fbaB - fructose-bisphosphate aldolase
    'G_b3115': 16.8,  # tpiA - triose-phosphate isomerase
    'G_b3236': 12.5,  # mdh - malate dehydrogenase
    'G_b3386': 9.2,   # rpe - ribulose-phosphate 3-epimerase
    'G_b3387': 8.8,   # rpiA - ribose-5-phosphate isomerase
    'G_b3603': 5.8,   # ppc - phosphoenolpyruvate carboxylase
    'G_b3731': 10.5,  # pgi - glucose-6-phosphate isomerase
    'G_b3732': 2.8,   # pfkB - phosphofructokinase
    'G_b3733': 11.2,  # gpmA - phosphoglycerate mutase
    'G_b3734': 13.5,  # eno - enolase
    'G_b3735': 8.9,   # pykF - pyruvate kinase
    'G_b3919': 3.2,   # glk - glucokinase
    'G_b3956': 7.5,   # pgl - 6-phosphogluconolactonase
    'G_b4015': 6.8,   # aceA - isocitrate lyase
    'G_b4016': 6.5,   # aceB - malate synthase A
    'G_b4025': 9.8,   # pgi - glucose-6-phosphate isomerase
    'G_b4090': 4.5,   # pfkA - phosphofructokinase
    'G_b4151': 8.2,   # maeB - malic enzyme
    'G_b4152': 11.5,  # talA - transaldolase
    'G_b4153': 12.8,  # tktB - transketolase
    'G_b4154': 13.2,  # zwf - glucose-6-phosphate dehydrogenase
    'G_b4232': 10.8,  # glk - glucokinase
    'G_b4395': 12.5,  # frdD - fumarate reductase
}

print(f"Gene expression data for {len(gene_expression)} genes")

# 1. Standard FBA (baseline)
print("\n1. Standard FBA (no expression):")
fba_sol = FBA(model)
print(f"   Growth rate: {fba_sol.fobj:.4f}")

# 2. GIMME (minimizes inconsistency with expression)
print("\n2. GIMME (Gene Inactivity Moderated by Metabolism and Expression):")
gimme_sol = GIMME(model, gene_expression, cutoff=25, growth_frac=0.9)
print(f"   Growth rate: {gimme_sol.fobj:.4f}")
print(f"   Active reactions: {sum(1 for v in gimme_sol.values.values() if abs(v) > 1e-6)}")

# 3. E-Flux (constrains bounds by expression)
print("\n3. E-Flux (Expression-based Flux):")
eflux_sol = eFlux(model, gene_expression)
print(f"   Growth rate: {eflux_sol.fobj:.4f}")
print(f"   Active reactions: {sum(1 for v in eflux_sol.values.values() if abs(v) > 1e-6)}")

# Compare key fluxes
print("\n4. Comparison of key fluxes:")
print(f"{'Reaction':<20} {'FBA':>10} {'GIMME':>10} {'E-Flux':>10}")
print("-" * 55)

key_rxns = ['R_PGI', 'R_PFK', 'R_PYK', 'R_CS', 'R_ACONTa', 'R_CYTBD']
for rxn in key_rxns:
    if rxn in fba_sol.values:
        print(f"{rxn:<20} {fba_sol.values.get(rxn, 0):>10.4f} "
              f"{gimme_sol.values.get(rxn, 0):>10.4f} "
              f"{eflux_sol.values.get(rxn, 0):>10.4f}")

print("\nDone! Expression-integrated FBA completed successfully.")
