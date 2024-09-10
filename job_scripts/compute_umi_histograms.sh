#!/bin/bash

#SBATCH --job-name=umi_histograms
#SBATCH --partition=all
#SBATCH --ntasks=15

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder>"
    echo "[-d <docker_image_path>] [-s <script_folder>]"
    exit 1
}

job_scripts_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$job_scripts_directory")"

# Variables to hold arguments
input_folder=""
output_folder=""

docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
script_folder="$repository_path"/scripts

# Parse command line arguments
while getopts ":i:o:d:s:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        s )
            script_folder=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] ||
[ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run UMI counting in docker
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "python3 /script_folder/compute_umi_histogram.py \
--input_folder /input_folder --output_folder /output_folder; \
chmod 777 -R /output_folder"
