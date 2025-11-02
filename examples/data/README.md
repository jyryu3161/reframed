# Gene Essentiality Data Files

이 디렉토리에는 gene knockout 실험 데이터가 포함되어 있습니다.

## 파일 설명

### 1. `ecoli_gene_essentiality.csv`
E. coli K-12 MG1655 유전자 essentiality 데이터 (CSV 형식)

**데이터 소스:**
- Keio collection (Baba et al., 2006, Mol Syst Biol)
- PEC database (Gerdes et al., 2003, J Bacteriol)

**성장 조건:** M9 minimal medium + glucose (aerobic, 37°C)

**컬럼:**
- `gene_id`: Gene identifier (모델의 gene ID와 일치)
- `essentiality`: essential, non-essential, conditionally-essential
- `confidence`: high, medium, low (실험 증거의 신뢰도)
- `condition`: Growth condition
- `reference`: 데이터 출처

**포맷 예시:**
```csv
gene_id,essentiality,confidence,condition,reference
G_b0008,non-essential,high,minimal_glucose,Baba2006
G_b0114,essential,high,minimal_glucose,Baba2006
G_b0875,essential,high,minimal_glucose,Baba2006
```

### 2. `ecoli_gene_essentiality.json`
E. coli K-12 MG1655 유전자 essentiality 데이터 (JSON 형식)

CSV와 동일한 데이터에 추가 정보 포함:
- 유전자 이름 및 설명
- Growth rate ratio (knockout/wild-type)
- 대사 경로 정보
- 상세한 metadata

**포맷 예시:**
```json
{
  "metadata": {
    "organism": "Escherichia coli K-12 MG1655",
    "condition": "minimal_glucose_aerobic",
    "references": ["Baba et al. (2006)"]
  },
  "genes": {
    "G_b0875": {
      "name": "pykF",
      "product": "pyruvate kinase I",
      "essentiality": "essential",
      "growth_rate_ratio": 0.04,
      "pathway": "glycolysis"
    }
  }
}
```

## 사용 방법

### Python에서 로딩

```python
from gene_essentiality_evaluation import load_essentiality_data

# CSV 로딩
ess_data = load_essentiality_data('data/ecoli_gene_essentiality.csv')

# JSON 로딩 (더 많은 정보 포함)
ess_data = load_essentiality_data('data/ecoli_gene_essentiality.json')

# 특정 조건만 필터링 (CSV만)
ess_data = load_essentiality_data(
    'data/ecoli_gene_essentiality.csv',
    condition='minimal_glucose'
)
```

### 통계

**E. coli core model (77 genes):**
- Essential genes: 24 (31.2%)
- Non-essential genes: 53 (68.8%)

**경로별 essentiality:**
- Glycolysis: 100% essential (11/11 genes)
- TCA cycle: 100% essential (1/1 genes)
- Pentose phosphate pathway: 80% essential (4/5 genes)
- Pyruvate metabolism: 100% essential (2/2 genes)

## Essentiality 기준

**Essential (필수):**
- Knockout 시 growth rate < 5% of wild-type
- 세포 생존에 필수적인 유전자

**Non-Essential (비필수):**
- Knockout 시 growth rate ≥ 80% of wild-type
- 제거해도 정상 성장 가능

**Conditionally-Essential (조건 의존):**
- 5% ≤ growth rate < 80% of wild-type
- 특정 조건에서만 필수

## 커스텀 데이터 추가

### 1. CSV 포맷

```csv
gene_id,essentiality,confidence,condition,reference
G_your_gene,essential,high,your_condition,YourStudy2025
```

**필수 컬럼:**
- `gene_id`: 모델의 gene ID와 일치해야 함
- `essentiality`: essential, non-essential, conditionally-essential 중 하나
- `condition`: 성장 조건 이름
- `confidence`: high, medium, low 중 하나

### 2. JSON 포맷

```json
{
  "metadata": {
    "organism": "Your organism",
    "condition": "your_condition",
    "description": "Description of your data"
  },
  "genes": {
    "G_your_gene": {
      "name": "gene_name",
      "essentiality": "essential",
      "confidence": "high",
      "growth_rate_ratio": 0.02
    }
  }
}
```

## 참고 문헌

1. **Baba T, et al. (2006)**
   Construction of Escherichia coli K-12 in-frame, single-gene knockout mutants: the Keio collection.
   *Molecular Systems Biology* 2:2006.0008

2. **Gerdes SY, et al. (2003)**
   Experimental determination and system level analysis of essential genes in Escherichia coli MG1655.
   *Journal of Bacteriology* 185:5673-5684

3. **Joyce AR, et al. (2006)**
   Experimental and computational assessment of conditionally essential genes in Escherichia coli.
   *Journal of Bacteriology* 188:8259-8271

## 다른 organism 데이터

동일한 포맷으로 다른 organism의 essentiality 데이터를 만들 수 있습니다:

- **Bacillus subtilis**: Kobayashi et al. (2003) PNAS
- **Mycobacterium tuberculosis**: Sassetti et al. (2003) PNAS
- **Saccharomyces cerevisiae**: Giaever et al. (2002) Nature
- **Pseudomonas aeruginosa**: Jacobs et al. (2003) PNAS

## 라이센스

이 예시 데이터는 published research를 기반으로 합니다.
실제 연구에 사용할 경우 원본 논문을 인용해주세요.
