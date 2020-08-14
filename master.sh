cd /project/6007297/vivek22/FM_RNA/FASTQ
# fix file names:
for file in *; do mv "$file" "$(echo "$file" | cut -c60-)"; done
for file in *; do mv -v "$file" "${file/.unmapped/}"; done
for file in *; do mv -v "$file" "${file/.fa/.Fa}"; done
cat -> fastqc.sh << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=400M
#SBATCH --time=00:40:00
module load fastqc/0.11.9
fastqc -t 10 *.fastq.gz
EOF
sbatch fastqc.sh
mkdir fastqc
mv *fastqc* fastqc
mv *slurm* fastqc
cd fastqc
multiqc .
mv /project/6007297/vivek22/FM_RNA/FASTQ/fastqc/ /project/6007297/vivek22/FM_RNA/
mkdir /project/6007297/vivek22/FM_RNA/FASTQ/z
mv /project/6007297/vivek22/FM_RNA/FASTQ/A00266_0043_1_i* /project/6007297/vivek22/FM_RNA/FASTQ/z
cd ..
mkdir salmon
cd salmon
# Download the reference transcriptome from the Gencode's FTP server
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
# prepare metadata (https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/):
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.v34.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
# make index file
cat -> salmon_index.sh << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10:59:00
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.3.0
salmon index -t gentrome.fa.gz -d decoys.txt -p 40 -i salmon_index --gencode
EOF
sbatch salmon_index.sh
################# Salmon quant script for all FASTQs (pwd = salmon) #################
for fq in /project/6007297/vivek22/FM_RNA/FASTQ/*.Fastq.gz; do
base=`basename $fq .fq`
cat -> ${base}.sh << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:59:00
# run salmon
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.3.0
salmon quant -i salmon_index             \
 -l A                                    \
 -r ${fq}                                \
 -p 40                                   \
 -o ${base}.salmon                       \
 --seqBias                               \
 --gcBias                                \
 --posBias                               \
 --numBootstraps 30
EOF
done
for sh in *Fastq.gz.sh ; do
sbatch ${sh}
done
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
