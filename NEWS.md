# futurervpa 0.0.0.9000

## rfuture & Fref

- 2012. 7. 28.; calc.rel.abund導入。SPR、YPRを計算する関数。
- 2012. 8. 1; 加入の関数を大幅改善
- 0.3f; Frec1, Frec2, Frecの計算オプションをfuture.vpaに追加
- 0.4: スケトウ対応
- 0.5: 2013.1.18 全種対応のためのバージョン
- 0.5c: 対馬のサバに対応するため、将来予測の部分にtimestepを導入
        資源量の単位はトン、尾数の単位は百万尾
- 2014.7.4; 黒田さん指摘により，Popeの引数を追加
- 2016.6.12; vpaのほうとファイルを分ける→future1.0
- MSY計算関数(SR.est)を追加
- 2017.08.10 SR.estにssbとRだけを単純に与えて再生産関係をフィットする機能を追加。引数SSB.datとR.datを与えると、他のデータは無視して,与えたデータにフィットさせる
- 2017.08.28 将来予測関数(future.vpa)にBbanオプションを追加
- 2017.11.24 ref.Fにmin.ageを追加
- 2017.11.24 ref.F, future.vpaにwaa.catch引数（漁獲量計算時のwaa。通常のwaa引数は資源量・親魚資源量計算時に使われる）。また、vpaへの引数にcaa.catchを特別に与えている場合には、将来予測や管理基準値、MSYの漁獲量の計算時にvpaに与えたcaa.catchを使うように変更
- 2017.12.25 Frecの引数に "scenario"と"Frange"オプションを追加
- scenario="catch.mean"とすると、将来のstochastic simulationにおける平均漁獲量Blimitで指定した値と一致するようになる
- Frange=c(0.01,2) などと指定する。探索するFの範囲。特に、scenario="catch.mean"を使う場合、解が２つ出てくるので、
- Frange=c(0.01,Fmsy), c(Fmsy,2) のように、２種類のFの範囲で探索することになる
- 2018.6.5 user interfaceの大幅改善。SR.estをいくつかの要素に分ける
- 2018.6.8, 11:05;  plot.kobe, plot.kobemat, get.kobematなど、結果のプロット系の関数を追加
- 2018.6.8, 16:30;  SR.fitの出力に加入の予測値が入るように変更
- 2018.6.11, get.kobemat2, plot.kobemat2を追加
- 2018.7.12 2.0にバージョンアップ。future.vpaで使っていない関数などを整理し、ベクトル化。
  加入関数で動くのはHS.rec(自己相関なし。残差のランダムリサンプリング対応), HS.recAR, BH.recAR, RI.recAI（←自己相関に対応。残差のランダムリサンプリングには未対応）のみ
- 2018.9.18 2.1; 岡村さんのABC計算用プログラムに対応するためのfuture.vpaに置き換え
                 資源量の面プロットを作る関数plotBfishを追加
                 岡村さん作成のest.MSY2（ARの考慮あり）、calc.abc関数を追加
- 2018.10.17; github上にアップ
