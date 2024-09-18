#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -g <gse>"
    exit 1
}

pipeline_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$pipeline_directory")"

# Variables to hold arguments
gse=""

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
slurm_log_folder="$repository_path"/slurm_logs
aligned_folder="/data/public/jkoubele/sc_data/aligned"
split_by_cell_folder="/data/public/jkoubele/sc_data/split_by_cell"


# Parse command line arguments
while getopts ":g:" opt; do
    case ${opt} in
        g )
            gse=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$gse" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

for sub_folder in "$aligned_folder"/"$gse"/*; do
  sample_name=$(basename "$sub_folder")
  echo "Submitting sample $sample_name"
  output_folder="$split_by_cell_folder"/"$gse"/"$sample_name"
  mkdir "$output_folder" -p
  sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
  "$repository_path"/job_scripts/split_bam_by_barcode.sh \
  -i "$sub_folder" -o "$output_folder" -r "$repository_path"
done
