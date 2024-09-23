#!/bin/bash

#SBATCH --job-name=call_errors
#SBATCH --partition=all
#SBATCH --ntasks=3

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -l <output_folder_locations>"
    echo "[-r <repository_path>] [-g <genome_folder] [-f <fasta_file_name>]"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
repository_path="/data/public/jkoubele/transcription-mistakes"
genome_folder="/data/public/jkoubele/reference_genomes/GRCm39"
fasta_file_name="Mus_musculus.GRCm39.dna.primary_assembly.fa"


# Parse command line arguments
while getopts ":i:o:r:g:f:l:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        r )
            repository_path=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        f )
            fasta_file_name=$OPTARG
            ;;
        l )
            output_folder_locations=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$repository_path" ] ||
[ -z "$genome_folder" ] || [ -z "$fasta_file_name" ] || [ -z "$output_folder_locations" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi


docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
script_folder="$repository_path"/scripts

docker load -i "$docker_image_path"
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
-v "$output_folder_locations":/output_folder_locations \
-v "$script_folder":/script_folder \
-v "$genome_folder":/genome_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "python3 /script_folder/call_errors.py \
--input_folder /input_folder --output_folder /output_folder --output_folder_locations /output_folder_locations \
--reference_genome_fasta_file /genome_folder/$fasta_file_name; \
chmod 777 -R /output_folder; chmod 777 -R /output_folder_locations"
