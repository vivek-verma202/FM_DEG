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
# Download the reference transcriptome from the FTP server
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz
# make index file
cat -> salmon_index.sh << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:59:00
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.1.0
salmon index \
-t ./Homo_sapiens.GRCh38.cdna.all.fa \
-p 10 \
-i salmon_index \
-k 31
EOF
sbatch salmon_index.sh
salloc --time=6:59:0 --ntasks=40 --mem-per-cpu=1G --account=def-ldiatc
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.1.0
for fq in /project/6007297/vivek22/FM_RNA/FASTQ/*.Fastq.gz
do
# create a prefix
base=`basename $fq .fq`
# run salmon
salmon quant -i salmon_index \
 -l A \
 -r $fq \
 -p 40 \
 -o $base.salmon \
 --seqBias \
 --useVBOpt \
 --numBootstraps 30 \
 --gcBias \
 --validateMappings

done
EOF
sbatch salmon_quant.sh
# check if files remaining (from FASTQ dir): 
diff  <(ls -1 /project/6007297/vivek22/FM_RNA/FASTQ/ | sed s/.Fastq.gz//g) <( ls -1 ./ | sed s/.Fastq.gz.salmon//g) 

salloc --time=2:59:0 --ntasks=10 --mem-per-cpu=500M --account=def-ldiatc
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.1.0
for fq in /project/6007297/vivek22/FM_RNA/FASTQ/rem/*.Fastq.gz
do
# create a prefix
base=`basename $fq .fq`
# run salmon
salmon quant -i salmon_index \
 -l A \
 -r $fq \
 -p 10 \
 -o $base.salmon \
 --seqBias \
 --useVBOpt \
 --numBootstraps 30 \
 --gcBias \
 --validateMappings
done
