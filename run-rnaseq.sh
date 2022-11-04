#Download the reference from NCBI genome (https://www.ncbi.nlm.nih.gov/genome/)
# 1) whole sequence (.fna file)
# 2) genome info (.gff file)


# Unzip raw data
gunzip ./raw_data/*

# Trim adapter or low quality sequences in local laptop using Terminal followed by Quality check by Fastqc
# https://github.com/FelixKrueger/TrimGalore
mkdir -p ./raw_data/trimmed
trim_galore -o ./raw_data/trimmed ./raw_data/*.fastq --paired --fastqc --illumina
mkdir -p ./raw_data/trimmed/seq_quality

# Save the quality info
mv ./raw_data/trimmed/*.txt ./raw_data/trimmed/seq_quality
mv ./raw_data/trimmed/*.html ./raw_data/trimmed/seq_quality
mv ./raw_data/trimmed/*.zip ./raw_data/trimmed/seq_quality

# Make an index genome by bowtie2
# https://github.com/BenLangmead/bowtie2
for file in ./*.fna; do cp $file $file.fa; done
#cp *.fna *.fa
mkdir -p ./raw_data/trimmed/Index
bowtie2-build *.fa ./raw_data/trimmed/Index/index


### Alignment of the sequence reads to the index genome by bowtie2
for f1 in ./raw_data/trimmed/*R1_val_1.fq ; do for f2 in ${f1%%R1_val_1.fq}"R2_val_2.fq" ; do bowtie2 -x ./raw_data/trimmed/Index/index -1 $f1 -2 $f2 -S ./$f1.sam ; done; done


### File format change by samtools
# https://github.com/samtools/samtools
for file in ./raw_data/trimmed/*.sam; do samtools view -bS $file -o ./$file.bam ; done
for file in ./raw_data/trimmed/*.bam; do samtools sort $file -o ./$file.sort.bam ; done

### Read counts by htseq
# https://htseq.readthedocs.io/en/release_0.11.1/index.html
mkdir -p ./raw_data/trimmed/count
for file in ./raw_data/trimmed/*.sort.bam ; do htseq-count -t gene -i Name --stranded=no -f bam -r pos ./$file ./*.gff > ./$file.count ; done
mv ./raw_data/trimmed/*.count ./raw_data/trimmed/count


