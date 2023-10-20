
mkdir data
cd data

# リファレンス配列の準備
## リファレンスのデータをダウンロード
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar xvzf Homo_sapiens_UCSC_hg38.tar.gz

mkdir hg38_noalt
cd hg38_noalt

# chr1-22, X, Y, Mにする
wget https://kero.hgc.jp/book/DSTEP_ProgramII_input_2020/data/script/reference.pl

perl reference.pl ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa genome.fa

grep "^>" genome.fa # 確認用

# ファイルの実行ディレクトリについてはきちんと整理を行う


# データの階層構造の作成
cd ~/RNAseq_analyzer/
mkdir tools
mkdir analysis_rnaseq
cd data
mkdir fastq


# そのほかのツールのインストール
cd ../tools
mkdir bin
export PATH=~/RNAseq_analyzer/tools/bin:$PATH


## FASTQC(ver. 0.12.1)
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod 750 FastQC/fastqc
ln -s ~/RNAseq_analyzer/tools/FastQC/fastqc bin
fastqc -h


## trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
sudo apt install openjdk-13-jre-headless
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar

# samtools
wget https://sourceforge.net/projects/samtools/files/samtools/1.18/samtools-1.18.tar.bz2
tar -jxvf samtools-1.18.tar.bz2
cd samtools-1.18/
./configure --prefix=/home/yukato/RNAseq_analyzer/tools  
make
# make clean
sudo make install
./samtools
make cd ..


wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip
unzip bowtie2-2.4.2-sra-linux-x86_64.zip
ln -s ~/RNAseq_analyzer/tools/bowtie2-2.4.2-sra-linux-x86_64/bowtie2* bin
bowtie2 -h


## STAR (mappingのツール)
wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
tar -xzf 2.7.11a.tar.gz
ln ~/RNAseq_analyzer/tools/STAR-2.7.11a/bin/Linux_x86_64/* bin
STAR

## featureCounts (カウントツール)
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
tar -xvzf subread-2.0.6-Linux-x86_64.tar.gz
ln -s ~/RNAseq_analyzer/tools/subread-2.0.6-Linux-x86_64/bin/featureCounts
  bin
featureCounts

## RSEM
conda install -c bioconda RSEM

## RSEMの準備
rsem-prepare-reference \
--num-threads 16 \
--gtf  ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf \
~/RNAseq_analyzer/data/hg38_noalt/genome.fa  \
~/data/RSEM_Reference


## Bowtieの準備
cd ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/
export PATH=~/RNAseq_analyzer/tools/bin:$PATH
bowtie2-build -f humRibosomal.fa humRibosomal

## STARの準備
cd ~/RNAseq_analyzer/data/hg38_noalt
mkdir STAR_index
cd STAR_index
ln -s ../genome.fa 
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles genome.fa --sjdbGTFfile ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.g
tf --sjdbOverhang 100



