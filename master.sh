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
salloc --time=2:59:0 --ntasks=40 --mem-per-cpu=2G --account=def-ldiatc
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.3.0
salmon index -t gentrome.fa.gz -d decoys.txt -p 40 -i salmon_index --gencode

################# test one (pwd = salmon) #############################
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.3.0
salmon quant -i salmon_index                           \
 -l A                                                  \
 -r /project/6007297/vivek22/FM_RNA/FASTQ/E05.Fastq.gz \
 -p 40                                                 \
 -o test                                               \
 --seqBias                                             \
 --gcBias                                              \
 --posBias                                             \
 --numBootstraps 30
 
################# script for all (pwd = salmon) #############################
for fq in /project/6007297/vivek22/FM_RNA/FASTQ/*.Fastq.gz; do
cat -> salmon_quant_${fq} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:59:00
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.3.0
salmon quant -i gencode.v34_salmon_1.3.0 \
 -l A                                    \
 -r ${fq}                                \
 -p 40                                   \
 -o ${fq}.salmon                         \
 --seqBias                               \
 --gcBias                                \
 --posBias                               \
 --numBootstraps 30
EOF

########## check scripts ##########
lt 
cat salmon_quant_...

########## submit scripts ##########
for fq in /project/6007297/vivek22/FM_RNA/FASTQ/*.Fastq.gz; do
sbatch salmon_quant_${fq}
done




>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>old_below
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
