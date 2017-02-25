# 近傍力のみを考慮した粒子運動学計算のCUDA実装
近傍力としてLJポテンシャルによる力を採用．


##MD-CPU[完成]
CPUでの実装，テスト・比較用

##MD-GPU-Naive[完成]
GPUでの実装，ナイーブ．移しただけ．

##MD-GPU-GRID[完成]
カット半径によるGRID化．GRID化はナイーブ．雛形．

##MD-GPU-GRID-SL
MD-GPU-GRIDをもとにカスタマイズ中


#可視化
![GAZO](https://i.gyazo.com/1a092037a4e5298ee195dc603cac48b1.png)


