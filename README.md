#how2use_try  

入力データの作成例./sparse_nonsymm の引数は,seed,10000行,1000列,各行20要素(seedに依存する乱数),出力名.  
gcc -o sparse_nonsymm sparse_nonsymm.c  
./sparse_nonsymm 89 10000 1000 20 sparse_data.txt  
make clean  
make  
\#下の引数は、得る特異対, s(sparse) or d(dense),accuracy = 10^(-x)  
./fileinput_gkl 30 s 7 < sparse_data.txt

\#SystemDの時は  
make -f Makefile_4_SystemD  
\#する。  
  
0927追記  
現在はfileinput_gkl.cが外部ファイル読み込みに非対応(ベンチマーク).使い方は  
make  
./fileinput_gkl 30 s 100 < dd1000000_100000_100  
  
fileinput_gkl.c.cが外部ファイル読み込み対応版  

# TRLAN
