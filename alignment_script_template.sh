#!/bin/bash

#SBATCH --job-name=align
#SBATCH --partition=all
#SBATCH --ntasks=15



# Variables to hold arguments
sra_folder="/data/public/jkoubele/SRA/FASTQ_SRR/"
output_folder="/data/public/jkoubele/sc_data/aligned/GSE/GSM/"
repo_folder="/data/public/jkoubele/transcription-mistakes/"
genome_folder="/data/public/jkoubele/reference_genomes/GRCm39/"
docker_image_path="/data/public/jkoubele/pol-II-analysis/docker_images/bioinfo_tools.tar"


# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run STAR aligner
docker run --rm \
-v "$sra_folder":/sra_folder \
-v "$output_folder":/output_folder \
-v "$repo_folder":/repo_folder \
-v "$genome_folder":/genome_folder \
--security-opt seccomp=unconfined bioinfo_tools /bin/sh -c \
"STAR \
--soloType CB_UMI_Simple \
--genomeDir /genome_folder/STAR_index \
--soloCBwhitelist /repo_folder/barcodes/BARCODES_FILE \
--soloUMIlen 12 \
--readFilesIn READS_2 \
READS_1 \
--readFilesCommand zcat \
--outFileNamePrefix /output_folder/ \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 50000000000 \
--outSAMattributes CB UB GX GN \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloFeatures GeneFull; \
chmod 777 -R /output_folder"
