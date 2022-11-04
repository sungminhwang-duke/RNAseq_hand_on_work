# Download the reference from NCBI genome (https://www.ncbi.nlm.nih.gov/genome/)
#   and place the files below in the "same directory where raw sequencing files"!! 
# 1) whole sequence (.fna file)
# 2) genome info (.gff file)


# Unzip raw data
gunzip ./Raw_materials/*

# Trim adapter or low quality sequences in local laptop using Terminal followed by Quality check by Fastqc
# https://github.com/FelixKrueger/TrimGalore
mkdir -p ./Raw_materials/trimmed
trim_galore -o ./Raw_materials/trimmed ./Raw_materials/*.fastq --paired --fastqc --illumina

# Save the quality info
mkdir -p ./Raw_materials/trimmed/seq_quality
mv ./Raw_materials/trimmed/*.txt ./Raw_materials/trimmed/seq_quality
mv ./Raw_materials/trimmed/*.html ./Raw_materials/trimmed/seq_quality
mv ./Raw_materials/trimmed/*.zip ./Raw_materials/trimmed/seq_quality

# Make an index genome by bowtie2
# https://github.com/BenLangmead/bowtie2
for file in ./Raw_materials/*.fna; do cp $file $file.fa; done
#cp *.fna *.fa
mkdir -p ./Raw_materials/trimmed/Index
bowtie2-build ./Raw_materials/*.fa ./Raw_materials/trimmed/Index/index


### Alignment of the sequence reads to the index genome by bowtie2
for f1 in ./Raw_materials/trimmed/*R1_val_1.fq ; do for f2 in ${f1%%R1_val_1.fq}"R2_val_2.fq" ; do bowtie2 -x ./Raw_materials/trimmed/Index/index -1 $f1 -2 $f2 -S ./$f1.sam ; done; done


### File format change by samtools
# https://github.com/samtools/samtools
for file in ./Raw_materials/trimmed/*.sam; do samtools view -bS $file -o ./$file.bam ; done
for file in ./Raw_materials/trimmed/*.bam; do samtools sort $file -o ./$file.sort.bam ; done

### Read counts by htseq
# https://htseq.readthedocs.io/en/release_0.11.1/index.html
mkdir -p ./Raw_materials/trimmed/count
for file in ./Raw_materials/trimmed/*.sort.bam ; do htseq-count -t gene -i Name --stranded=no -f bam -r pos ./$file ./Raw_materials/*.gff > ./$file.count ; done
mv ./Raw_materials/trimmed/*.count ./Raw_materials/trimmed/count


