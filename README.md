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

### ディレクトリを作成する
```
mkdir RNAseq
mkdir Reference
```

### インストール色々
mambaが入っている前提で行う
```
mamba install -c bioconda fastqc -y
mamba install -c bioconda multiqc -y
mamba install -c bioconda star -y
mamba install -c bioconda rsem -y
mamba install -c bioconda samtools -y
```

### reference配列とgtfファイルをインストールと準備
- インストールもとは、Ensembl [https://asia.ensembl.org/Homo_sapiens/Info/Index] 
- バージョンが変わるので、定期的に見直すと良いかも(執筆時2024年2月2日)
```
cd Refence
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.111.gtf.gz
```

### STARのindexの作成
- 解析のルートディレクトリで実行
```
STAR --runMode genomeGenerate \
--genomeDir ~/Reference/STAR_index \
--runThreadN 16 \
--genomeFastaFiles ~/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/Reference/Homo_sapiens.GRCh38.111.gtf
```

### RSEMのindexの作成
- 解析のルートディレクトリで実行
```
rsem-prepare-reference \
--num-threads 16 \
--gtf ~/Reference/Homo_sapiens.GRCh38.111.gtf \
~/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
~/Reference/RSEM_Reference
```

### 解析の実行
RNAseqフォルダに解析対象のfastqファイルを入れて以下をコマンドラインから実行する
```
cd RNAseq # 適宜変更
bash pipeline.sh
```





# 参考資料
- 東大式 生命データサイエンス即戦力講座〜ゲノム、エピゲノム、トランスクリプトームからシングルセルまで、大規模データ解析で論文を書くためのR&Pythonツールボックス
- windowsでRNA-seq解析 [https://zenn.dev/rchiji/books/cd3bc4612d79b4]

