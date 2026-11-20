# spatialfamily — 空間統計と空間計量経済学（jamovi モジュール）

都市経済学・地域経済学・空間経済学で使う空間分析：座標、境界ファイル、隣接リストから空間ウェイト行列を作り、空間的集積を検定・地図化し、空間回帰モデルを最尤法で推定します。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：Spatial）

| 分析 | 内容 |
|---|---|
| **空間ウェイト** | W を作成して記述：リンク数、1 単位あたりの隣接数、孤立単位、連結成分、接続図、隣接数のヒストグラム。入力：座標列（k 近傍、距離帯、逆距離；平面または経緯度）、ポリゴンの shapefile または GeoJSON（クイーン／ルーク隣接、または重心ベース）、隣接リスト（エッジリスト csv/txt、GAL、GWT）。高次隣接、対称化、行基準化／二値／全体基準化。 |
| **大域的自己相関（Moran、Geary、G）** | 複数の変数について Moran の I、Geary の C、Getis-Ord の大域 G；ランダム化・正規性・置換による推測；Moran 散布図。 |
| **局所的自己相関（LISA、Gi\*）** | 局所 Moran の I（High-High、Low-Low、High-Low、Low-High）、Getis-Ord Gi と Gi\*；条件付き置換 p 値、FDR または Bonferroni 補正；クラスター地図（ポリゴンまたは点）、Moran 散布図、単位ごとの表；クラス・統計量・p 値をデータセットに保存。 |
| **空間パネル回帰** | 均衡パネル：個体・時間・双方向の固定効果（またはプール）での SAR、SEM、SDM、SDEM、SLX の最尤推定（within 変換、T ln\|I - ρW\|）、直接・間接・総効果。 |
| **地理的加重回帰（GWR）** | 空間的に変わる係数：bi-square、ガウス、指数、tricube、箱型カーネル；適応型・固定バンド幅を AICc または交差検証で選択；局所標準誤差と t 値、局所 R²、係数の地図、OLS との AICc・有効パラメータ数の比較。 |
| **産業集積** | 地域×産業の雇用表から：立地ジニ係数、ハーフィンダール指数、Ellison-Glaeser γ（事業所数を指定可）、地域の Krugman 特化指数；立地係数をデータセットに保存。 |
| **重力モデル（PPML）** | ゼロを含む二国間フローのポアソン擬似最尤、高次元固定効果（出発地、目的地、ペア、出発地×期、目的地×期、期）；頑健・クラスター・二元クラスター標準誤差。 |
| **空間 DID（波及）** | 双方向固定効果 DID、処置単位の周りの距離リング、空間ラグ処置、共変量；単位でクラスター、頑健、または Conley 空間 HAC 標準誤差。 |
| **空間回帰** | OLS、SLX、SAR（ラグ）、SEM（誤差）、SDM（ダービン）、SDEM、SAC を最尤法で推定；OLS 残差の LM 検定（誤差、ラグ、頑健、SARMA）；平均の直接・間接・総効果（シミュレーション標準誤差）；残差の Moran の I；7 モデルを対数尤度・AIC・BIC で比較；予測値と残差の保存。 |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

## ウェイトの入力元
- **座標**：数値列 2 つ。経緯度の場合は「経度／緯度」にチェック（大圏距離、km）。
- **境界ファイル**：ポリゴンの shapefile（`.shp` と同じフォルダの `.dbf`）または GeoJSON。パスを入力し、データの id 列と、ファイル内の id フィールド名を指定します。
- **隣接リストファイル**：1 行 1 リンク（`from,to` または `from,to,weight`）、または GAL／GWT ファイル。id はデータの id 列と一致させてください。
- モジュール付属のサンプルは `sample:columbus.geojson`、`sample:columbus_edges.csv`、`sample:nc.shp`、`sample:nc.geojson` で指定できます。

## 検証
エンジンは base R と Matrix のみで書かれ、R パッケージ `spdep`、`spatialreg` と照合済みです（`tests/reference.R`）：Moran の I、Geary の C、大域 G の分散と p 値、局所 Moran（Ii、E、Var、z）、Gi と Gi\*、LM 検定、SAR/SEM/SDM/SDEM/SAC/SLX の最尤推定（対数尤度、係数、rho、lambda、観測情報標準誤差）と直接・間接効果が表示精度で一致します。`tests/weights_ref.R` は shapefile／GeoJSON の読み込み、クイーン／ルーク隣接、高次隣接、k 近傍、距離帯、GAL の読み込みを `spdep` と `sf` に照合します。さらに `splm`、`GWmodel`、`fixest` および Stata 19（`spregress`、`spxtregress`、`ppmlhdfe`、`reghdfe`）と一致することを確認済みです。サンプルデータ：Columbus 市の近隣犯罪（49 単位）、ノースカロライナ州の郡別 SIDS（100 郡）、およびシミュレーションの空間パネル、地域別産業雇用、二国間貿易、波及 DID。

## 制限
空間パネルは均衡パネルと固定効果のみ（変量効果・動学パネルなし）。GWR は 4000 単位まで、非定常性のモンテカルロ検定はありません。空間 RD は `rdfamily` の地理 RD を使ってください。Duranton-Overman の距離ベース集積指標（事業所座標が必要）は含まれません。1500 単位を超える空間回帰は行基準化ウェイトが必要です。jamovi にはインタラクティブな地図がないため地図は静止画です。ポリゴンの shapefile と GeoJSON のみ読み込めます。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *spatialfamily* (jamovi module). https://github.com/mingchao-sun/spatialfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/spatialfamily/issues
