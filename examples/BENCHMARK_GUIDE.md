# 벤치마크 프레임워크 가이드 (Benchmark Framework Guide)

## 🎯 빠른 시작 (Quick Start)

### 1. 기본 벤치마크 실행

가장 간단한 방법으로 여러 방법론을 비교하려면:

```bash
cd examples
python quick_benchmark.py
```

**결과:**
- 터미널에 요약 결과 출력
- `benchmark_results.csv`: 메트릭 요약
- `benchmark_detailed_fluxes.csv`: 모든 반응의 flux 값

---

## 📊 결과 해석 (Understanding Results)

### 벤치마크 요약 테이블

```
Method                    Status    Growth  Active Rxns  Total Flux  Time (s)
------------------------------------------------------------------------------
FBA                       OPTIMAL   0.8739          48      518.42     0.008
GIMME(cutoff=25,gf=0.9)   OPTIMAL  62.6783          52      472.05     0.013
E-Flux                    OPTIMAL   0.0097          52        6.11     0.007
```

**컬럼 설명:**
- **Method**: 방법론 이름 및 매개변수
- **Status**: 최적화 상태 (OPTIMAL이 정상)
- **Growth**: 성장률 (biomass flux)
- **Active Rxns**: 0이 아닌 flux를 가진 반응 수
- **Total Flux**: 모든 flux의 절대값 합 (대사 활동 수준)
- **Time**: 실행 시간 (초)

### 주요 지표 (Key Metrics)

1. **Growth Rate** (성장률)
   - 높을수록 좋음 (일반적으로)
   - FBA는 최대 성장률 제공
   - GIMME는 growth_frac에 따라 조정됨

2. **Active Reactions** (활성 반응 수)
   - 사용된 대사 경로의 복잡도
   - 적을수록 단순한 대사 경로
   - 발현 기반 방법은 경로를 제한함

3. **Total Flux** (총 flux)
   - 대사 활동의 전반적 수준
   - 높을수록 더 많은 대사 활동

4. **Execution Time** (실행 시간)
   - 계산 효율성
   - 대규모 분석 시 중요

---

## 🔧 커스텀 벤치마크 작성

### 예제 1: 특정 방법론만 비교

```python
from reframed import load_cbmodel
from benchmark_transcriptomics import (
    TranscriptomicsBenchmark,
    GIMMEMethod,
    generate_sample_expression
)

# 모델 로딩
model = load_cbmodel('tests/data/e_coli_core.xml.gz')
gene_exp = generate_sample_expression(model, 'aerobic_glucose')

# 벤치마크 생성
benchmark = TranscriptomicsBenchmark(model, gene_exp)

# GIMME 방법론만 다양한 매개변수로 비교
for cutoff in [10, 25, 50, 75, 90]:
    benchmark.add_method(GIMMEMethod(cutoff=cutoff, growth_frac=0.9))

# 실행
results = benchmark.run_all()
benchmark.print_summary()

# 결과 내보내기
benchmark.export_results('gimme_sensitivity.csv')
```

### 예제 2: 여러 조건에서 비교

```python
conditions = ['aerobic_glucose', 'anaerobic_glucose', 'aerobic_acetate']
all_results = {}

for condition in conditions:
    gene_exp = generate_sample_expression(model, condition)
    benchmark = TranscriptomicsBenchmark(model, gene_exp)

    benchmark.add_methods([
        FBAMethod(),
        GIMMEMethod(cutoff=25, growth_frac=0.9),
        EFluxMethod(),
    ])

    results = benchmark.run_all(verbose=False)
    all_results[condition] = results

    # 조건별 결과 저장
    benchmark.export_results(f'results_{condition}.csv')
```

### 예제 3: Flux 비교 분석

```python
# 특정 반응들의 flux 비교
key_reactions = [
    'R_GLCpts',  # Glucose uptake
    'R_PGI',     # Glycolysis
    'R_PFK',     # Glycolysis
    'R_CS',      # TCA cycle
    'R_CYTBD',   # Respiration
]

benchmark.compare_fluxes(key_reactions)

# 또는 가장 변동이 큰 반응들 자동 선택
benchmark.compare_fluxes([], top_n=20)
```

---

## 🆕 신규 방법론 추가 (Adding New Methods)

### Step 1: 템플릿 선택

`custom_method_template.py`에는 4가지 템플릿이 있습니다:

1. **SimpleExpressionMethod**: 발현 기반 제약 추가
2. **OptimizationBasedMethod**: 커스텀 목적 함수
3. **TwoPhaseMethod**: 이단계 최적화
4. **MLBasedMethod**: 머신러닝 통합

### Step 2: 메서드 구현

```python
from benchmark_transcriptomics import TranscriptomicsMethod, BenchmarkResult
import time

class MyCustomMethod(TranscriptomicsMethod):
    """내 새로운 방법론"""

    def __init__(self, param1=1.0, param2=0.5):
        super().__init__(f"MyMethod(p1={param1},p2={param2})")
        self.param1 = param1
        self.param2 = param2

    def run(self, model, gene_exp, **kwargs):
        try:
            start_time = time.time()

            # 여기에 알고리즘 구현
            # 1. 발현 데이터 처리
            # 2. 모델 수정 또는 최적화 문제 정의
            # 3. FBA 또는 커스텀 최적화 실행

            # 예: 간단한 FBA
            from reframed import FBA
            solution = FBA(model)

            exec_time = time.time() - start_time

            return self._calculate_metrics(
                solution,
                self.name,
                exec_time,
                {'param1': self.param1, 'param2': self.param2}
            )

        except Exception as e:
            return BenchmarkResult(
                method_name=self.name,
                status="ERROR",
                growth_rate=0.0,
                execution_time=0.0,
                active_reactions=0,
                total_flux=0.0,
                error=str(e)
            )
```

### Step 3: 벤치마크에 추가

```python
benchmark = TranscriptomicsBenchmark(model, gene_exp)

# 기존 방법론과 함께 비교
benchmark.add_methods([
    FBAMethod(),
    GIMMEMethod(cutoff=25, growth_frac=0.9),
    MyCustomMethod(param1=2.0, param2=0.8),  # 내 방법론!
])

results = benchmark.run_all()
benchmark.print_summary()
```

---

## 📈 매개변수 최적화 (Parameter Optimization)

### Grid Search 방식

```python
# 두 개의 매개변수를 동시에 탐색
cutoffs = [10, 25, 50, 75]
growth_fracs = [0.7, 0.8, 0.9, 0.95]

benchmark = TranscriptomicsBenchmark(model, gene_exp)
benchmark.add_method(FBAMethod())  # Baseline

for cutoff in cutoffs:
    for gf in growth_fracs:
        benchmark.add_method(
            GIMMEMethod(cutoff=cutoff, growth_frac=gf)
        )

results = benchmark.run_all(verbose=False)
benchmark.export_results('parameter_grid_search.csv')

# CSV를 pandas로 분석
import pandas as pd
df = pd.DataFrame([
    {
        'cutoff': r.parameters.get('cutoff'),
        'growth_frac': r.parameters.get('growth_frac'),
        'growth_rate': r.growth_rate,
        'active_reactions': r.active_reactions,
    }
    for r in results if 'GIMME' in r.method_name
])

# 최적 매개변수 찾기
best_params = df.loc[df['growth_rate'].idxmax()]
print(f"Best parameters: cutoff={best_params['cutoff']}, "
      f"growth_frac={best_params['growth_frac']}")
```

---

## 🔍 고급 분석 (Advanced Analysis)

### 1. 방법론 간 상관관계 분석

```python
import numpy as np
from scipy.stats import pearsonr

# flux 데이터 추출
methods = ['FBA', 'GIMME(cutoff=25,gf=0.9)', 'E-Flux']
flux_data = {}

for result in results:
    if result.method_name in methods:
        flux_data[result.method_name] = result.flux_values

# 상관계수 계산
reactions = list(flux_data['FBA'].keys())
for m1 in methods:
    for m2 in methods:
        if m1 < m2:
            fluxes1 = [flux_data[m1].get(r, 0) for r in reactions]
            fluxes2 = [flux_data[m2].get(r, 0) for r in reactions]

            corr, pval = pearsonr(fluxes1, fluxes2)
            print(f"{m1} vs {m2}: r={corr:.3f} (p={pval:.3e})")
```

### 2. 대사 경로별 분석

```python
# 주요 경로별로 반응 그룹화
pathways = {
    'Glycolysis': ['R_PGI', 'R_PFK', 'R_FBA', 'R_TPI', 'R_GAPD', 'R_PGK', 'R_PGM', 'R_ENO', 'R_PYK'],
    'TCA': ['R_CS', 'R_ACONTa', 'R_ACONTb', 'R_ICDHyr', 'R_AKGDH', 'R_SUCOAS', 'R_SUCD1i', 'R_FUM', 'R_MDH'],
    'Respiration': ['R_CYTBD', 'R_NADH16', 'R_ATPS4r'],
}

for pathway_name, reactions in pathways.items():
    print(f"\n=== {pathway_name} ===")
    benchmark.compare_fluxes(reactions)
```

### 3. 반응 활성화/비활성화 패턴

```python
# 각 방법론에서 활성화된 반응 집합
for result in results:
    if result.status == 'OPTIMAL':
        active = {r for r, v in result.flux_values.items() if abs(v) > 1e-6}
        print(f"{result.method_name}: {len(active)} active reactions")

# 모든 방법론에서 활성화된 core 반응
all_active = [
    {r for r, v in result.flux_values.items() if abs(v) > 1e-6}
    for result in results if result.status == 'OPTIMAL'
]

core_reactions = set.intersection(*all_active)
print(f"\nCore active reactions (all methods): {len(core_reactions)}")
```

---

## 💡 팁과 모범 사례 (Tips & Best Practices)

### 1. 발현 데이터 품질 확인

```python
import numpy as np

# 기본 통계
values = list(gene_exp.values())
print(f"Expression stats:")
print(f"  Mean: {np.mean(values):.2f}")
print(f"  Std:  {np.std(values):.2f}")
print(f"  Min:  {np.min(values):.2f}")
print(f"  Max:  {np.max(values):.2f}")
print(f"  Coverage: {len(gene_exp)}/{len(model.genes)} genes")

# 분포 확인
import matplotlib.pyplot as plt
plt.hist(values, bins=50)
plt.xlabel('Expression Level')
plt.ylabel('Frequency')
plt.title('Gene Expression Distribution')
plt.savefig('expression_distribution.png')
```

### 2. 계산 성능 최적화

- 큰 모델은 SCIP 대신 Gurobi/CPLEX 사용
- 매개변수 grid search 시 병렬 처리 고려
- `verbose=False`로 출력 최소화

### 3. 결과 검증

```python
# 1. Growth rate가 합리적인지 확인
assert 0 < result.growth_rate < 10, "Unusual growth rate"

# 2. Mass balance 확인
from reframed import FVA
fva_result = FVA(model, obj_frac=0.99)

# 3. 실험 데이터와 비교
experimental_fluxes = {
    'R_GLCpts': 10.0,
    'R_PYK': 2.5,
    # ...
}

for rxn, exp_flux in experimental_fluxes.items():
    pred_flux = result.flux_values.get(rxn, 0)
    error = abs(pred_flux - exp_flux) / exp_flux
    print(f"{rxn}: predicted={pred_flux:.2f}, "
          f"experimental={exp_flux:.2f}, error={error*100:.1f}%")
```

### 4. 재현 가능성

```python
# Random seed 고정
import numpy as np
np.random.seed(42)

# 발현 데이터 저장
import json
with open('gene_expression.json', 'w') as f:
    json.dump(gene_exp, f, indent=2)

# 벤치마크 설정 저장
config = {
    'model': 'e_coli_core.xml.gz',
    'condition': 'aerobic_glucose',
    'methods': [
        {'name': 'GIMME', 'cutoff': 25, 'growth_frac': 0.9},
        {'name': 'E-Flux'},
    ],
    'date': '2025-11-02',
}

with open('benchmark_config.json', 'w') as f:
    json.dump(config, f, indent=2)
```

---

## ❓ 자주 묻는 질문 (FAQ)

**Q: 왜 GIMME의 growth rate가 FBA보다 높나요?**

A: GIMME의 목적 함수는 성장률 최대화가 아니라 발현 불일치 최소화입니다. 출력된 growth_rate는 실제로 목적 함수 값일 수 있습니다. biomass flux를 확인하려면:

```python
biomass_flux = result.flux_values[model.biomass_reaction]
```

**Q: 여러 조건의 결과를 하나의 테이블로 정리하려면?**

```python
import pandas as pd

all_data = []
for condition in conditions:
    gene_exp = generate_sample_expression(model, condition)
    benchmark = TranscriptomicsBenchmark(model, gene_exp)
    # ... run benchmark ...

    for result in benchmark.results:
        all_data.append({
            'Condition': condition,
            'Method': result.method_name,
            'Growth': result.growth_rate,
            'Active_Rxns': result.active_reactions,
        })

df = pd.DataFrame(all_data)
df.to_csv('multi_condition_results.csv', index=False)
```

**Q: 새 방법론이 ERROR 상태로 나옵니다.**

A: 에러 메시지를 확인하세요:

```python
for result in results:
    if result.error:
        print(f"{result.method_name}: {result.error}")
```

일반적인 원인:
- 발현 데이터 형식 오류
- 모델 제약 조건 충돌
- 잘못된 매개변수 값

---

## 📚 추가 자료 (Additional Resources)

- **기본 예제**: `simple_gimme_eflux.py`
- **종합 분석**: `ecoli_gimme_eflux_example.py`
- **빠른 벤치마크**: `quick_benchmark.py`
- **전체 프레임워크**: `benchmark_transcriptomics.py`
- **커스텀 템플릿**: `custom_method_template.py`

### 논문 참고

- GIMME: Becker & Palsson (2008) PLoS Comp Biol
- E-Flux: Colijn et al. (2009) PLoS Comp Biol
- FBA: Orth et al. (2010) Nature Biotech
- ReFramed: https://reframed.readthedocs.io

---

## 🤝 기여 및 피드백 (Contributing)

새로운 방법론을 개발하셨나요? 벤치마크 결과를 공유하고 싶으신가요?

1. 커스텀 방법론을 `custom_method_template.py` 스타일로 작성
2. 벤치마크 결과를 CSV로 내보내기
3. GitHub Issues에 공유

Happy benchmarking! 🚀
