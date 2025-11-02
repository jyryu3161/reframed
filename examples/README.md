# ReFramed Examples

이 디렉토리에는 ReFramed 라이브러리를 사용한 다양한 대사 모델링 예제가 포함되어 있습니다.

This directory contains various examples demonstrating the use of ReFramed for metabolic modeling.

## 📁 파일 목록 (Files)

### 1. `simple_gimme_eflux.py` - 간단한 시작 예제
**난이도:** ⭐ 초급 (Beginner)

전사체 데이터를 대사 모델과 통합하는 가장 간단한 예제입니다.

- **포함 내용:**
  - 대장균 코어 모델 로딩
  - 샘플 유전자 발현 데이터
  - GIMME 분석 실행
  - E-Flux 분석 실행
  - 결과 비교

- **실행 방법:**
  ```bash
  cd examples
  python simple_gimme_eflux.py
  ```

### 2. `quick_benchmark.py` - 빠른 벤치마크 🆕
**난이도:** ⭐⭐ 중급 (Intermediate)

여러 방법론을 동시에 비교하는 벤치마크 예제입니다.

- **포함 내용:**
  - FBA, pFBA, GIMME, E-Flux 동시 실행
  - 자동 성능 메트릭 계산
  - Flux 비교 및 시각화
  - 결과 CSV 내보내기

- **실행 방법:**
  ```bash
  cd examples
  python quick_benchmark.py
  ```

- **결과 파일:**
  - `benchmark_results.csv` - 요약 통계
  - `benchmark_detailed_fluxes.csv` - 상세 flux 데이터

### 3. `benchmark_transcriptomics.py` - 종합 벤치마크 프레임워크 🆕
**난이도:** ⭐⭐⭐ 고급 (Advanced)

확장 가능한 벤치마크 프레임워크로 여러 예제를 포함합니다.

- **포함 내용:**
  - 기본 방법론 비교
  - 매개변수 민감도 분석
  - 다중 조건 분석
  - 상세 메트릭 및 통계
  - 결과 내보내기 기능

- **실행 방법:**
  ```bash
  cd examples
  python benchmark_transcriptomics.py
  ```

- **주요 클래스:**
  - `TranscriptomicsBenchmark` - 벤치마크 실행 클래스
  - `TranscriptomicsMethod` - 방법론 기본 클래스
  - `BenchmarkResult` - 결과 저장 클래스

### 4. `custom_method_template.py` - 신규 방법론 템플릿 🆕
**난이도:** ⭐⭐⭐ 고급 (Advanced)

새로운 전사체 통합 방법론을 구현하기 위한 템플릿 코드입니다.

- **포함된 템플릿:**
  1. **SimpleExpressionMethod** - 간단한 발현 기반 방법
  2. **OptimizationBasedMethod** - 최적화 기반 방법
  3. **TwoPhaseMethod** - 이단계 최적화 방법
  4. **MLBasedMethod** - 머신러닝 기반 방법 (구조)

- **사용 방법:**
  ```python
  # 1. 템플릿 복사 및 수정
  # 2. run() 메서드 구현
  # 3. 벤치마크에 추가
  from custom_method_template import SimpleExpressionMethod

  benchmark.add_method(SimpleExpressionMethod(threshold=2.0))
  ```

### 5. `gene_essentiality_evaluation.py` - Gene Essentiality 평가 모듈 🆕
**난이도:** ⭐⭐⭐ 고급 (Advanced)

Gene knockout 실험 데이터로 모델을 검증하는 모듈입니다.

- **주요 기능:**
  - Essentiality 데이터 로딩 (CSV/JSON)
  - Gene knockout 시뮬레이션
  - 예측 vs 실험 데이터 비교
  - 성능 메트릭 계산 (Accuracy, Precision, Recall, F1, MCC)
  - Confusion matrix 분석

- **사용 방법:**
  ```python
  from gene_essentiality_evaluation import (
      load_essentiality_data,
      GeneEssentialityEvaluator
  )

  # 데이터 로딩
  ess_data = load_essentiality_data('data/ecoli_gene_essentiality.csv')

  # 평가
  evaluator = GeneEssentialityEvaluator(model, ess_data)
  evaluation = evaluator.evaluate_fba()

  # 결과 확인
  evaluator.print_metrics(evaluation)
  ```

### 6. `essentiality_benchmark_example.py` - Essentiality 벤치마크 예제 🆕
**난이도:** ⭐⭐ 중급 (Intermediate)

Gene essentiality 데이터로 모델을 평가하는 완전한 예제입니다.

- **포함 내용:**
  - FBA/pFBA로 knockout 시뮬레이션
  - 실험 데이터와 비교
  - Confusion matrix 및 메트릭
  - Threshold 민감도 분석
  - 경로별 essentiality 분석

- **실행 방법:**
  ```bash
  cd examples
  python essentiality_benchmark_example.py
  ```

- **결과 파일:**
  - `essentiality_predictions.csv` - 각 유전자별 예측
  - `essentiality_comparison.csv` - 방법론 비교

### 7. `ecoli_gimme_eflux_example.py` - 상세 분석 예제
**난이도:** ⭐⭐⭐ 고급 (Advanced)

GIMME와 E-Flux 방법론에 대한 포괄적인 분석 예제입니다.

- **포함 내용:**
  - 조건 특이적 유전자 발현 시뮬레이션
  - 호기성/혐기성 조건 비교
  - 다양한 탄소원에 대한 분석
  - 매개변수 민감도 분석
  - 상세한 결과 시각화

- **실행 방법:**
  ```bash
  cd examples
  python ecoli_gimme_eflux_example.py
  ```

## 🧬 GIMME vs E-Flux 방법론 비교

### GIMME (Gene Inactivity Moderated by Metabolism and Expression)

**원리:**
- 발현량이 낮은 유전자와 연관된 반응에 페널티 부여
- 최소 성장률을 유지하면서 페널티 최소화
- 이단계 최적화: (1) 최소 성장 달성, (2) 낮은 발현 반응 최소화

**장점:**
- 생물학적으로 타당한 flux 분포 생성
- 성장률 제약 조건 명시적 설정 가능
- 발현 데이터의 노이즈에 상대적으로 강건함

**적용 사례:**
- 조직 특이적 대사 모델링
- 암세포 대사 분석
- 조건 특이적 대사 경로 예측

**참고문헌:**
Becker, S. A., & Palsson, B. Ø. (2008). Context-specific metabolic networks are consistent with experiments. *PLoS Computational Biology*, 4(5), e1000082.

### E-Flux (Expression-based Flux)

**원리:**
- 유전자 발현 수준에 비례하여 반응의 상한 제약
- 발현량을 최대값으로 정규화
- 제약된 모델로 FBA 수행

**장점:**
- 구현이 간단하고 직관적
- 계산 속도가 빠름
- 발현 데이터를 직접적으로 반영

**적용 사례:**
- 전사체 데이터 통합
- 대규모 스크리닝 연구
- 빠른 조건별 대사 예측

**참고문헌:**
Colijn, C., Brandes, A., Zucker, J., Lun, D. S., Weiner, B., Farhat, M. R., ... & Galagan, J. E. (2009). Interpreting expression data with metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid production. *PLoS Computational Biology*, 5(8), e1000489.

## 🔬 벤치마크 프레임워크 사용법 (Benchmark Framework Usage) 🆕

### 기본 벤치마크 실행

```python
from reframed import load_cbmodel
from benchmark_transcriptomics import (
    TranscriptomicsBenchmark,
    FBAMethod,
    GIMMEMethod,
    EFluxMethod,
    generate_sample_expression
)

# 1. 모델과 발현 데이터 준비
model = load_cbmodel('ecoli_model.xml')
gene_exp = generate_sample_expression(model, 'aerobic_glucose')

# 2. 벤치마크 생성
benchmark = TranscriptomicsBenchmark(model, gene_exp)

# 3. 비교할 방법론 추가
benchmark.add_methods([
    FBAMethod(),
    GIMMEMethod(cutoff=25, growth_frac=0.9),
    EFluxMethod(),
])

# 4. 실행
results = benchmark.run_all(verbose=True)

# 5. 결과 확인
benchmark.print_summary()
benchmark.compare_fluxes(['R_PGI', 'R_PFK', 'R_PYK'])

# 6. 결과 내보내기
benchmark.export_results('results.csv')
```

### 신규 방법론 추가하기

**Step 1: 템플릿 선택**

```python
from benchmark_transcriptomics import TranscriptomicsMethod, BenchmarkResult
```

**Step 2: 클래스 구현**

```python
class MyNewMethod(TranscriptomicsMethod):
    def __init__(self, parameter1=1.0):
        super().__init__(f"MyMethod(p1={parameter1})")
        self.parameter1 = parameter1

    def run(self, model, gene_exp, **kwargs):
        # 알고리즘 구현
        # ...
        return self._calculate_metrics(solution, self.name, exec_time, params)
```

**Step 3: 벤치마크에 추가**

```python
benchmark.add_method(MyNewMethod(parameter1=2.0))
```

### 매개변수 민감도 분석

```python
# GIMME cutoff 값 비교
for cutoff in [10, 25, 50, 75, 90]:
    benchmark.add_method(GIMMEMethod(cutoff=cutoff, growth_frac=0.9))

results = benchmark.run_all()
benchmark.print_summary()
```

### 다중 조건 분석

```python
conditions = ['aerobic_glucose', 'anaerobic_glucose', 'aerobic_acetate']

for condition in conditions:
    gene_exp = generate_sample_expression(model, condition)
    benchmark = TranscriptomicsBenchmark(model, gene_exp)

    benchmark.add_methods([
        FBAMethod(),
        GIMMEMethod(),
        EFluxMethod(),
    ])

    print(f"\n=== Condition: {condition} ===")
    benchmark.run_all()
    benchmark.print_summary()
```

## 📊 사용 워크플로우 (Workflow)

### 1. 데이터 준비
```python
# RNA-seq 또는 마이크로어레이 데이터를 딕셔너리로 변환
gene_expression = {
    'gene_id_1': expression_value_1,
    'gene_id_2': expression_value_2,
    # ...
}
```

### 2. 모델 로딩
```python
from reframed import load_cbmodel

model = load_cbmodel('path/to/model.xml')
```

### 3. GIMME 실행
```python
from reframed import GIMME

solution = GIMME(
    model,
    gene_expression,
    cutoff=25,        # 발현 임계값 백분위수 (25th percentile)
    growth_frac=0.9   # 최소 성장 비율 (야생형의 90%)
)
```

### 4. E-Flux 실행
```python
from reframed import eFlux

solution = eFlux(
    model,
    gene_expression,
    scale_rxn='R_GLCpts',  # 정규화 기준 반응 (선택사항)
    scale_value=10.0        # 정규화 값 (선택사항)
)
```

### 5. 결과 분석
```python
# 성장률 확인
print(f"Growth rate: {solution.fobj}")

# 주요 flux 확인
for rxn_id in ['R_PGI', 'R_PFK', 'R_PYK']:
    print(f"{rxn_id}: {solution.values[rxn_id]:.4f}")

# 활성 반응 개수
active_reactions = sum(1 for v in solution.values.values() if abs(v) > 1e-6)
print(f"Active reactions: {active_reactions}")
```

## 🔧 매개변수 가이드

### GIMME 매개변수

| 매개변수 | 기본값 | 설명 | 권장 범위 |
|---------|--------|------|-----------|
| `cutoff` | 25 | 발현 임계값 백분위수 | 10-50 |
| `growth_frac` | 0.9 | 최소 성장 비율 | 0.5-0.99 |
| `parsimonious` | False | 절약적 최적화 사용 | True/False |

**Tips:**
- `cutoff`이 낮을수록 더 많은 반응에 페널티
- `growth_frac`이 높을수록 FBA와 유사한 결과
- `parsimonious=True`는 flux 분포를 더 단순하게 만듦

### E-Flux 매개변수

| 매개변수 | 기본값 | 설명 |
|---------|--------|------|
| `scale_rxn` | None | 정규화 기준 반응 |
| `scale_value` | 1.0 | 정규화 값 |
| `parsimonious` | False | pFBA 사용 |

**Tips:**
- `scale_rxn`으로 실험적 uptake rate에 맞출 수 있음
- 발현 데이터는 0보다 큰 값이어야 함

## 📖 추가 자료

### 공식 문서
- [ReFramed Documentation](https://reframed.readthedocs.io)
- [API Reference](https://reframed.readthedocs.io/en/latest/api.html)

### 관련 논문
1. **GIMME**: Becker & Palsson (2008) - Context-specific metabolic networks
2. **E-Flux**: Colijn et al. (2009) - Interpreting expression data with metabolic flux models
3. **FBA**: Orth et al. (2010) - What is flux balance analysis?

### 대사 모델 데이터베이스
- [BiGG Models](http://bigg.ucsd.edu/) - 고품질 genome-scale 모델
- [BioModels](https://www.ebi.ac.uk/biomodels/) - SBML 형식 모델 저장소
- [MetaNetX](https://www.metanetx.org/) - 대사 네트워크 통합 플랫폼

## ❓ FAQ

**Q: 모델에 없는 유전자의 발현 데이터는 어떻게 되나요?**
A: 자동으로 무시됩니다. 모델에 있는 유전자만 사용됩니다.

**Q: 발현 데이터가 없는 유전자는 어떻게 처리되나요?**
A: GIMME와 E-Flux 모두 발현 데이터가 없는 유전자는 중간 수준으로 가정하거나 제약하지 않습니다.

**Q: GIMME와 E-Flux 중 무엇을 선택해야 하나요?**
A:
- 생물학적 타당성 중시 → GIMME
- 계산 속도 중시 → E-Flux
- 노이즈가 많은 데이터 → GIMME
- 대규모 스크리닝 → E-Flux

**Q: 음수 발현값이 있으면 어떻게 하나요?**
A: 발현값은 항상 0 이상이어야 합니다. 전처리 과정에서 정규화하세요.

**Q: 다른 대사 모델 형식도 사용 가능한가요?**
A: ReFramed는 SBML 형식을 지원합니다. COBRApy 모델도 변환 가능합니다.

## 🐛 문제 해결

### 일반적인 오류

**1. "Infeasible solution"**
```python
# 원인: 너무 제한적인 제약
# 해결: growth_frac 낮추기
solution = GIMME(model, gene_exp, growth_frac=0.7)  # 0.9 → 0.7
```

**2. "Gene not found in model"**
```python
# 원인: 유전자 ID 불일치
# 해결: 모델의 유전자 ID 확인
print(list(model.genes.keys())[:10])  # 첫 10개 유전자 ID 확인
```

**3. "Expression values must be positive"**
```python
# 원인: 음수 발현값
# 해결: 데이터 전처리
gene_exp = {k: max(0.1, v) for k, v in gene_exp.items()}
```

## 💡 도움말

문제가 발생하거나 질문이 있으면:
- GitHub Issues: https://github.com/cdanielmachado/reframed/issues
- Documentation: https://reframed.readthedocs.io

## 📝 라이센스

이 예제들은 ReFramed 라이브러리와 동일한 Apache 2.0 라이센스를 따릅니다.
