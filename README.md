# RNA-seqの解析

1. fastqからcount, TPMを算出する
2. その後の解析
  a. pydeseq2
  b. pygsea
  c. ComplexHeatmapによる可視化

# 1. fastqからcount, TPMを算出する

## 実行するライブラリ
- FastQC
- fastp
- STAR
- RSEM

## a. 実行環境の構築
- pythonのパッケージインストーラーのmambaを利用する
- 環境構築に関しては、各OSで異なる場合があるので一つずつ行うのが良い。

### step 1: ディレクトリを作成する
```
mkdir RNAseq
mkdir RNAseq/Reference
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

Mac miniで行う場合STAR(2.7.11a)はエラーが出てしまうので、STAR(ver2.7.5a)を使用するとエラー回避できる。
mambaが2.7.11b以上のバージョンをインストール可能になれば、そちらのバージョンを使用した方が良いかもしれない。

### step3: reference配列とgtfファイルをインストールと準備
- インストールもとは、Ensembl [https://asia.ensembl.org/Homo_sapiens/Info/Index] 
- バージョンが変わるので、定期的に見直すと良いかも(執筆時2024年2月2日)
```
cd Reference
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.111.gtf.gz
```

### step:4 STARのindexの作成
- 解析のルートディレクトリで実行
```
STAR --runMode genomeGenerate \
--genomeDir ~/RNAseq/Reference/STAR_index \
--runThreadN 16 \
--genomeFastaFiles ~/RNAseq/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/RNAseq/Reference/Homo_sapiens.GRCh38.111.gtf
```

### step5: RSEMのindexの作成
- 解析のルートディレクトリで実行
```
rsem-prepare-reference \
--num-threads 16 \
--gtf ~/RNAseq/Reference/Homo_sapiens.GRCh38.111.gtf \
~/RNAseq/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
~/RNAseq/Reference/RSEM_Reference
```

## b. 解析の実行
RNAseqフォルダに解析対象のfastqファイルを入れて以下をコマンドラインから実行する
(bash または、zshをまなんで自作してみてください)
```
cd RNAseq # 適宜変更
bash pipeline.sh
```


## c. 備忘録

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
--readFilesIn <<your_fastq_1>> <<your_fastq_2>>
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
<<your-bam-file>> \
~/RNAseq/Reference/RSEM_Reference \
f_out
```

# 2. その後の解析
## 発現が変動している遺伝子を統計的な解析により明らかにしたい場合
### DESeq2を利用する場合
```
python3 pydeg2.py # ファイルを追加中です
```
### サンプル数が少ない場合、fold-changeで発現変動遺伝子を抽出したい場合
```
python3 extract_deg.py　# ファイルを追加中です
```

## enrichment解析を行いたい場合
### GSEApyを実行する
普通のGO, GSEA, ssGSEA, GSVAなどできるので便利です。
使用例は以下の通りです。


そのほかのデータベースを利用した情報は、私のwebサイトに記載しましたのでご参照下さい。

# 参考資料
- 東大式 生命データサイエンス即戦力講座〜ゲノム、エピゲノム、トランスクリプトームからシングルセルまで、大規模データ解析で論文を書くためのR&Pythonツールボックス (羊土社）
- windowsでRNA-seq解析 [https://zenn.dev/rchiji/books/cd3bc4612d79b4]
- PyDESeq2 [https://pydeseq2.readthedocs.io/en/latest/]
- GSEAPY [https://gseapy.readthedocs.io/en/latest/introduction.html]
