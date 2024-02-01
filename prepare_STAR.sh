# installに関しては、x86_64じゃないとだめかも
# arch # CPUの確認用
mamba install -c bioconda star -y

mkdir ref
cd ref
STAR --runMode genomeGenerate \
--genomeDir ~/ref/STAR_index \
--runThreadN 16 \
--genomeFastaFiles ~/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/Reference/Homo_sapiens.GRCh38.111.gtf
