# spatialfamily — 空间统计与空间计量（jamovi 模块）

面向城市、区域与空间经济学研究的空间分析：由坐标、边界文件或邻接表构造空间权重矩阵，检验并绘制空间集聚，并用极大似然估计空间回归模型。 英文完整说明见 [README.md](README.md)。

## 分析（jamovi 菜单：Spatial）

| 分析 | 内容 |
|---|---|
| **空间权重** | 构造并描述 W：连接数、每个单位的邻居数、孤立单位、连通分量、连接图和邻居数直方图。来源：坐标列（k 近邻、距离带、反距离；平面或经纬度）、多边形 shapefile 或 GeoJSON（Queen/Rook 邻接，或基于质心）、邻接表（边列表 csv/txt、GAL、GWT）；高阶邻接、对称化、行标准化/二值/全局标准化。 |
| **全局自相关（Moran、Geary、G）** | 可同时检验多个变量的 Moran's I、Geary's C 和 Getis-Ord 全局 G；随机化、正态或置换推断；Moran 散点图。 |
| **局部自相关（LISA、Gi\*）** | 局部 Moran's I（高-高、低-低、高-低、低-高），Getis-Ord Gi 与 Gi\*；条件置换 p 值，FDR 或 Bonferroni 校正；集聚地图（多边形或点）、Moran 散点图、逐单位表；可把类别、统计量和 p 值保存回数据集。 |
| **空间面板回归** | 平衡面板：个体、时间或双向固定效应（或混合）下的 SAR、SEM、SDM、SDEM、SLX 极大似然估计（组内变换，T ln\|I - ρW\|），直接/间接/总效应。 |
| **地理加权回归（GWR）** | 系数随空间变化：双平方、高斯、指数、三次立方、箱型核；自适应或固定带宽，按 AICc 或交叉验证选择；局部标准误与 t 值、局部 R²、系数地图，AICc 与有效参数个数对比 OLS。 |
| **产业集聚** | 由地区×产业就业表计算：区位基尼系数、赫芬达尔指数、Ellison-Glaeser γ（可用企业数）、地区 Krugman 专业化指数，区位熵可保存回数据集。 |
| **引力模型（PPML）** | 含零流量的双边流量泊松伪极大似然，高维固定效应（来源地、目的地、配对、来源地×期、目的地×期、期）；稳健、聚类和双向聚类标准误。 |
| **空间 DID（溢出）** | 双向固定效应 DID，处理单位周围的距离环带、空间滞后处理变量、协变量；按单位聚类、稳健或 Conley 空间 HAC 标准误。 |
| **空间回归** | 极大似然估计 OLS、SLX、SAR（空间滞后）、SEM（空间误差）、SDM（杜宾）、SDEM 和 SAC；对 OLS 残差的 LM 检验（误差、滞后、稳健、SARMA）；平均直接、间接和总效应（模拟标准误）；残差的 Moran's I；七个模型按对数似然、AIC、BIC 并列比较；可保存拟合值与残差。 |

每个分析的输出部分都有"方法说明"和"数据与设置指南"。

## 权重的来源
- **坐标**：两列数值；经纬度需勾选"经度/纬度"（球面距离，单位 km）。
- **边界文件**：多边形 shapefile（`.shp` 与同目录的 `.dbf`）或 GeoJSON；填写路径，选择数据中的 id 列，并填写文件中对应的 id 字段名。
- **邻接表文件**：每行一条连接（`from,to` 或 `from,to,weight`）或 GAL/GWT 文件；id 需与数据中的 id 列一致。
- 模块自带的示例文件可写作 `sample:columbus.geojson`、`sample:columbus_edges.csv`、`sample:nc.shp`、`sample:nc.geojson`。

## 验证
引擎只用基础 R 和 Matrix，已与 R 包 `spdep`、`spatialreg` 对照（`tests/reference.R`）：Moran's I、Geary's C、全局 G 及其方差与 p 值，局部 Moran（Ii、E、Var、z）、Gi 与 Gi\*，LM 检验，SAR/SEM/SDM/SDEM/SAC/SLX 的极大似然估计（对数似然、系数、rho、lambda、观测信息标准误）以及直接/间接效应，均在显示精度内一致。`tests/weights_ref.R` 将 shapefile 与 GeoJSON 读取、Queen/Rook 邻接、高阶邻接、k 近邻、距离带和 GAL 读取与 `spdep`、`sf` 对照。另与 `splm`、`GWmodel`、`fixest` 以及 Stata 19（`spregress`、`spxtregress`、`ppmlhdfe`、`reghdfe`）对照一致（`tests/`、`../stata_validation/07_spatial.do`）。示例数据：哥伦布市邻里犯罪（49 个单位）、北卡罗来纳县 SIDS（100 个县），以及模拟的空间面板、地区产业就业、双边贸易和溢出 DID 数据。

## 局限
空间面板需要平衡面板和固定效应（无随机效应、动态面板）；GWR 最多 4000 个单位，无非平稳性的蒙特卡洛检验；空间断点回归请用 `rdfamily` 的地理 RD；Duranton-Overman 距离型集聚指标（需要企业坐标）未包含。超过 1500 个单位时空间回归要求行标准化权重。jamovi 无交互式地图，地图为静态图；仅读取多边形 shapefile 与 GeoJSON。

## 安装
- **jamovi 模块库**：在 jamovi 右上角打开模块按钮，进入 *jamovi library*，搜索并安装（模块上架后）。
- **侧载**：下载 `.jmo` 文件，在同一菜单选 *Install from file*。`.jmo` 只适用于构建它的操作系统、处理器类型和 jamovi 版本。
- **从源码**：需要 R 和 `jmvtools`：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（路径按本机 jamovi 修改）。

## 引用与许可
引用：Sun, M. *spatialfamily* (jamovi module). https://github.com/mingchao-sun/spatialfamily 。请同时引用所用方法的原始文献（输出中列出）。

许可：见 `DESCRIPTION`（GPL）；第三方数据见 `NOTICE`。问题与建议：https://github.com/mingchao-sun/spatialfamily/issues
