# 設計書

## 概要

置換基付き共役コア構造探索システムは、CSDデータベースから置換基を取得し、共役コアに結合させて2D周期境界条件下で構造最適化を行い、STEP1の参照構造とのRMSDを評価するパイプラインシステムです。各計算ステップを自動化し、複数の置換基候補を効率的に処理できる設計とします。

## アーキテクチャ

システムは以下のモジュール構成とします：

```
substituted-conjugated-core-search/
├── data_acquisition/          # CSDデータ取得・スクリーニング
├── structure_builder/         # 構造結合・初期最適化
├── quantum_calculation/       # Gaussian・antechamber処理
├── periodic_optimization/     # 2D周期境界条件最適化
├── analysis/                  # RMSD計算・結果評価
├── workflow_manager/          # パイプライン制御
└── utils/                     # 共通ユーティリティ
```

## コンポーネントとインターフェース

### 1. データ取得モジュール (data_acquisition)

**責任:** CSDデータベースから置換基構造を取得・フィルタリング

**主要クラス:**
- `CSDConnector`: CSD Python APIを使用したデータベース接続・検索
- `SubstituentScreener`: 置換基のスクリーニング・フィルタリング
- `StructureExtractor`: 構造データの抽出・変換

**インターフェース:**
```python
class SubstituentData:
    structure_id: str
    molecular_formula: str
    coordinates: List[Atom]
    connection_points: List[int]

def get_substituents(screening_criteria: Dict) -> List[SubstituentData]
```

### 2. 構造構築モジュール (structure_builder)

**責任:** 共役コアと置換基の結合、初期構造最適化

**主要クラス:**
- `ConjugatedCore`: 共役コア構造の管理
- `BondingSiteSelector`: 結合位置の選択・管理
- `StructureCombiner`: 構造の結合処理
- `InitialOptimizer`: 簡単な構造最適化

**インターフェース:**
```python
class CombinedStructure:
    core_structure: ConjugatedCore
    substituents: List[SubstituentData]
    bonding_sites: List[int]
    optimized_coordinates: List[Atom]

def combine_structures(core: ConjugatedCore,
                      substituent: SubstituentData,
                      bonding_site: int) -> CombinedStructure
```

### 3. 量子計算モジュール (quantum_calculation)

**責任:** Gaussian計算とRESP電荷計算の自動化

**主要クラス:**
- `GaussianRunner`: Gaussian計算の実行・管理
- `AntechamberProcessor`: antechamber処理の実行
- `RESPChargeCalculator`: RESP電荷計算の統合処理

**インターフェース:**
```python
class QuantumResult:
    optimized_structure: CombinedStructure
    resp_charges: List[float]
    calculation_log: str

def calculate_resp_charges(structure: CombinedStructure) -> QuantumResult
```

### 4. 周期最適化モジュール (periodic_optimization)

**責任:** 2D周期境界条件下での構造最適化（STEP1コード再利用）

**主要クラス:**
- `PeriodicBoundarySetup`: 2D周期境界条件の設定
- `Step1CodeAdapter`: STEP1コードのアダプター
- `PeriodicOptimizer`: 周期境界条件下最適化の実行

**インターフェース:**
```python
class PeriodicResult:
    optimized_structure: CombinedStructure
    energy: float
    convergence_info: Dict

def optimize_with_periodic_boundary(structure: CombinedStructure,
                                   charges: List[float]) -> PeriodicResult
```

### 5. 解析モジュール (analysis)

**責任:** RMSD計算と結果評価

**主要クラス:**
- `RMSDCalculator`: 共役コア部分のRMSD計算
- `StructureComparator`: 構造比較・評価
- `ResultEvaluator`: 結果の評価・ランキング

**インターフェース:**
```python
class AnalysisResult:
    rmsd_value: float
    core_alignment: Dict
    evaluation_score: float

def calculate_core_rmsd(reference: ConjugatedCore,
                       optimized: PeriodicResult) -> AnalysisResult
```

### 6. ワークフロー管理モジュール (workflow_manager)

**責任:** パイプライン全体の制御・進捗管理

**主要クラス:**
- `PipelineController`: パイプライン全体の制御
- `ProgressTracker`: 進捗状況の追跡
- `ErrorHandler`: エラー処理・リカバリ

## データモデル

### 原子構造データ
```python
@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float
    charge: Optional[float] = None

@dataclass
class Bond:
    atom1_idx: int
    atom2_idx: int
    bond_order: float
```

### 計算設定データ
```python
@dataclass
class CalculationSettings:
    gaussian_method: str = "B3LYP"
    basis_set: str = "6-31G*"
    periodic_boundary_params: Dict
    optimization_tolerance: float = 1e-6
```

## エラーハンドリング

### エラー分類
1. **データ取得エラー**: CSD Python API接続失敗、データ不整合
2. **構造構築エラー**: 結合失敗、座標異常
3. **計算エラー**: Gaussian実行失敗、収束失敗
4. **ファイルI/Oエラー**: ファイル読み書き失敗

### エラー処理戦略
- 各モジュールで例外をキャッチし、適切なログを出力
- 計算失敗時は次の候補に進む（スキップ機能）
- 重要なエラーは処理を停止し、ユーザーに通知
- 中間結果を保存し、再開可能な設計

## テスト戦略

### 単体テスト
- 各モジュールの主要クラス・メソッドをテスト
- モックを使用した外部依存関係の分離
- 異常系のテストケースを含む

### 統合テスト
- モジュール間の連携をテスト
- 小規模なテストデータでパイプライン全体を実行
- STEP1コードとの統合テスト

### 受入テスト
- 実際の研究データを使用したエンドツーエンドテスト
- RMSD計算の精度検証
- パフォーマンステスト（複数候補の並列処理）

## 実装上の考慮事項

### STEP1コードの再利用
- 既存の2D周期境界条件コードをアダプターパターンで統合
- インターフェースの統一化
- 設定パラメータの互換性確保

### 計算効率化
- 置換基候補の並列処理
- 中間結果のキャッシュ機能
- 計算リソースの効率的な利用

### 拡張性
- 新しい置換基スクリーニング条件の追加
- 異なる量子計算ソフトウェアへの対応
- 結果評価指標の追加