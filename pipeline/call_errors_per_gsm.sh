#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -e <gse> -m <gsm> "
    exit 1
}

pipeline_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$pipeline_directory")"

# Variables to hold arguments
gse=""
gsm=""

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
slurm_log_folder="$repository_path"/slurm_logs
split_by_cell_folder="/data/public/jkoubele/sc_data/split_by_cell"
detected_errors_folder="/data/public/jkoubele/sc_data/detected_errors"


# Parse command line arguments
while getopts ":e:m:" opt; do
    case ${opt} in
        e )
            gse=$OPTARG
            ;;
        m )
            gsm=$OPTARG
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
if [ -z "$gse" ] || [ -z "$gsm" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

for sub_folder in "$split_by_cell_folder"/"$gse"/"$gsm"/*; do
  subset_name=$(basename "$sub_folder")
  echo "Submitting $subset_name"
  output_folder="$detected_errors_folder"/"$gse"/"$gsm"/"$subset_name"
  mkdir "$output_folder" -p
  sbatch --nodelist=beyer-n02 --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
  "$repository_path"/job_scripts/call_errors.sh \
  -i "$sub_folder" -o "$output_folder" -r "$repository_path"
done
