# fastqファイルがあるフォルダに移動する

STAR \
--genomeDir ~/ref/STAR_index \
--runThreadN 16 \
--outFileNamePrefix  your_file_name \
--quantMode TranscriptomeSAM \
--outSAMtype None \
--readFilesCommand zcat \
--readFilesIn your_file_1 your_file_2
