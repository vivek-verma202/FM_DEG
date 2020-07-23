cd /project/6007297/vivek22/FM_RNA/FASTQ
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
salloc --time=1:0:0 --mem=100G --account=def-ldiatc --verbose
module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.4  salmon/1.1.0
