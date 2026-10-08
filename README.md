# cp2k: 粘土エッジの原子電荷をCP2K(PBE、Hirshfeld)で作る

ピロフィライト(粘土)のエッジ付きリボンと、その上の水・Mg置換・Na⁺について、局所環境ごとの原子電荷を
CP2Kで計算する。SOAP → 電荷(のちに SOAP → χ → QEq)モデルの教師データが目的。

## 使い方(スパコン)
```bash
git clone https://github.com/haru2225/cp2k.git && cd cp2k
qsub run_cp2k.pbs
```
これだけ。`run_cp2k.pbs` が順に次を行う(終わった構造は `charges.dat` があればスキップ。再投入すれば続きから):
1. `runs/`(乾いたリボン7構造): 構造最適化 → 電荷計算。完了済みの5構造(`bulk`、`rib_x_o00_si`、`rib_y_o00_al`、`rib_y_o00_si`、`rib_y_o37_al`)は結果を同梱してあるのでスキップされ、
   残りの `rib_x_o00_al`(前回SCF不収束)と `rib_y_o37_si`(途中で停止)だけ計算する。
2. `runs_wet/`(リボン + 水 (+ Mg/Na)、12系): 与えた構造での1点計算のみ。

SCFが収束しない場合は自動でリトライする(最適化3回、1点計算2回)。1回目は速いOT、失敗したら直前の構造から、
対角化 + Broyden混合 + 300 K Fermi–Dirac smearing の頑健な設定で続ける。時間切れで止まった最適化も最後の構造から再開する。
smearingを使うと、エネルギーは僅かに変わる(電荷への影響は小さいはずだが未検証)。

`PBS` の指定: `-q sc16`、`select=1:ncpus=16:mpiprocs=16`、`walltime=24:00:00`(仮置き。キューの上限に合わせて `run_cp2k.pbs` の `#PBS` 行を編集)。
残りの見込み: 乾いた2構造(最適化は1〜2時間 + 電荷15分)、湿った12系(1系あたり十数分〜数十分)、合計で数時間〜10時間程度。

環境は自動で整える:
- センターのサンプル(`/home/center/app/CP2K/cp2k_20251.sh`)の `module load` / `source` / `export PATH|LD_LIBRARY_PATH...` の行を再生する(別ファイルは `qsub -v CP2K_ENV_SCRIPT=...`、無効化は `=none`)。
- CP2K実行ファイル(`cp2k.psmp` / `cp2k.popt`)と、`BASIS_MOLOPT` / `GTH_POTENTIALS` のあるデータフォルダを探す(`-v CP2K_EXE=...,CP2K_DATA_DIR=...` で指定可)。
- 計算ノードに python は不要(`tail` と `awk` だけ)。CP2Kの標準出力は `runs*/<名前>/*.stdout`、失敗時はジョブログに末尾を出す。
- 一部だけ計算: `qsub -v NAME=rib_x_o00_al run_cp2k.pbs`、片方だけ: `qsub -v RUNS=runs_wet run_cp2k.pbs`。

## 結果の取り出しと解析
```bash
tar czf results.tar.gz runs/*/{relaxed.xyz,sp.out,opt.out,charges.dat} runs_wet/*/{sp.out,charges.dat}
python3 analyze_charges.py     # 要 numpy, ase, dscribe, scikit-learn: SOAP電荷モデルの検証(学習に使わないエッジ向きで評価)
python3 analyze_wet.py         # 水(とMg/Na)によるリボン原子の電荷変化 dq = q_wet − q_dry (numpyのみ)
```

## 構造
- `runs/` 乾いたリボン(`structures/` から作った入力): `bulk`(3D周期)、y法線リボン4(`rib_y_o00_{si,al}`、`rib_y_o37_{si,al}`)、x法線リボン2(`rib_x_o00_{si,al}`)。
  `o00`/`o37` は切る位置、`si`/`al` はプロトンの配分(Si–OH優先/Al–OH₂優先)。形式電荷で中性。
- `runs_wet/` 湿った系(`structures_wet/`): `wet_<リボン>_w<k>` は水のみ(中性)、`wet_<リボン>_wNa<k>` は内部Al→Mg置換 + Na⁺ 1個(層電荷 −1、元のSi8Al3.5Mg0.5モデルと同密度)。k=1,2 は古典MDの6 ps、12 psのスナップショット。
  水は12〜23分子(約0.85 g/cm³)で、エッジ間の真空部に置いた。原子順は「リボン → Na → 水(O H H)」。
- 入力の条件: PBE、GTHポテンシャル、最適化は SZV-MOLOPT-SR-GTH / 300 Ry、電荷計算は DZVP-MOLOPT-SR-GTH / 400 Ry、全方向周期セル + 真空。

## スクリプト
| ファイル | 役割 |
|---|---|
| `run_cp2k.pbs` | 上記のジョブ(qsubはこれだけ) |
| `make_cp2k.py` | `structures*/` から `runs*/<名前>/` の入力を作り直す(生成済みをコミット済み) |
| `hirshfeld_to_charges.py` | `sp.out` → `charges.dat`(awk版と同等。ローカル用) |
| `build_edges.py` | `C2000.gro`(ClayCode、MIT)からリボンを作る |
| `build_wet.py` | 緩和したリボンに水・Mg・Na⁺を入れる(要 numpy と LAMMPS の python モジュール) |
| `analyze_charges.py`, `analyze_wet.py` | 解析 |

## 注意(確認済み/未確認)
- `run_cp2k.pbs` の流れは、偽のCP2K/mpirunによる試験と、ローカルのCP2K 2026.2での電荷計算(バルク、頑健設定)で確認した。**センターの実機のスパコンでは未実行**(`module` の動作や既定キューなど)。
- 頑健設定が構造最適化の不収束を実際に救うかは未確認。
- 湿った系のうち、Na/Mgを含む129原子の電荷計算をローカルで実行中(確認は未完了)。
- Hirshfeld電荷はClayFFの約1/3.8のスケール(バルク: Si +0.55、Al +0.42、O −0.27)。MDに使うにはスケール換算が必要。
- 水の配置は古典力場(ClayFF + SPC)の分布。Mgサイトは再緩和していない。
