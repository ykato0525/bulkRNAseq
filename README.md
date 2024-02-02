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

## 1. 実行環境の構築
- pythonのパッケージインストーラーのmambaを利用する
- 環境構築に関しては、各OSで異なる場合があるので一つずつ行うのが良い。

### step 1: ディレクトリを作成する
```
mkdir RNAseq
mkdir Reference
```

### step 2: インストール色々
mambaが入っている前提で行う
```
mamba install -c bioconda fastqc -y
mamba install -c bioconda multiqc -y
mamba install -c bioconda star -y
mamba install -c bioconda rsem -y
mamba install -c bioconda samtools -y
```

### step3: reference配列とgtfファイルをインストールと準備
- インストールもとは、Ensembl [https://asia.ensembl.org/Homo_sapiens/Info/Index] 
- バージョンが変わるので、定期的に見直すと良いかも(執筆時2024年2月2日)
```
cd Refence
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.111.gtf.gz
```

### step:4 STARのindexの作成
- 解析のルートディレクトリで実行
```
STAR --runMode genomeGenerate \
--genomeDir ~/Reference/STAR_index \
--runThreadN 16 \
--genomeFastaFiles ~/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/ref/Homo_sapiens.GRCh38.111.gtf
```

### step5: RSEMのindexの作成
- 解析のルートディレクトリで実行
```
rsem-prepare-reference \
--num-threads 16 \
--gtf ~/Reference/Homo_sapiens.GRCh38.111.gtf \
~/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
~/Reference/RSEM_Reference
```

## 2. 解析の実行
RNAseqフォルダに解析対象のfastqファイルを入れて以下をコマンドラインから実行する
```
cd RNAseq # 適宜変更
bash pipeline.sh
```

## 備忘録

### 個別のファイルの実行方法について
pipeline.shには、全てのファイルを一括で送信するためのコードを記載してあるが、ライブラリの個々のパラメータについて備忘録として記載する。

#### fastqc
```
fastqc -t 12 --nogroup <<fastqファイル>> ; fastqc -t 12 --nogroup <<fastqファイル>>
```
- 引数について
  -　-t : スレッド数
  -  --no-group: 数塩基単位でまとめて解析するのを防ぐ(基本的にいれた方が良い)
  -  -o: 保存されるディレクトリを指定する   

### STAR
```
STAR \
--genomeDir ~/Reference/STAR_index \
--runThreadN 16 \
--outFileNamePrefix SRR22571458_ \
--quantMode TranscriptomeSAM \
--outSAMtype None \
--readFilesCommand zcat \ # macの場合はgzcat
--readFilesIn fastq_1 fastq_2
```


### RSEM
```
rsem-calculate-expression \
-p 16 \
--paired-end \
--alignments \
--append-names \
--estimate-rspd \
--no-bam-output \
your"bam"file \
~/Reference/RSEM_Reference \
f_out
```



# 参考資料
- 東大式 生命データサイエンス即戦力講座〜ゲノム、エピゲノム、トランスクリプトームからシングルセルまで、大規模データ解析で論文を書くためのR&Pythonツールボックス
- windowsでRNA-seq解析 [https://zenn.dev/rchiji/books/cd3bc4612d79b4]

