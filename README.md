# bulk RNA-seqのカウント値の算出を自動化する

## 実行環境
WSL2上のUbuntu20.04を想定

## 実行するライブラリ
- FastQC
- multiQC
- fastp
- STAR
- RSEM
- biomart

## 実行環境の構築
pythonのパッケージインストーラーのmambaを利用する
### インストール色々
```
mamba install -c bioconda fastqc -y
mamba install -c bioconda multiqc -y
mamba install -c bioconda star -y
mamba install -c bioconda rsem -y
mamba install -c bioconda samtools -y
```
### reference配列とgtfファイルをインストールする
インストールもとは、Ensembl https://asia.ensembl.org/Homo_sapiens/Info/Index
ほかでも良いです。


### STARの準備




## 手順
パッケージのインストーラーのmambaを使用して、必要なツールのダウンロード


# 参考資料
東大式 生命データサイエンス即戦力講座〜ゲノム、エピゲノム、トランスクリプトームからシングルセルまで、大規模データ解析で論文を書くためのR&Pythonツールボックス
