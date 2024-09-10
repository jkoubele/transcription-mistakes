# Function to display help
show_help() {
  echo "Usage: $0 [-o output_folder] [--help]"
  echo ""
  echo "Options:"
  echo "  -o      Specify the output folder (default is 'repository_path/docker_images')."
  echo "  --help  Display this help message."
}

# Check for --help before parsing other options
for arg in "$@"; do
  case $arg in
    --help)
      show_help
      exit 0
      ;;
  esac
done

# default value of output_folder
script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
output_folder="$repository_path/docker_images"

# Parse command line options
while getopts ":o:" opt; do
  case ${opt} in
    o )
      output_folder=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-o output_folder]"
      exit 1
      ;;
  esac
done

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# build and save docker images
docker build -t bioinfo_tools "$repository_path/dockerfiles/bioinfo_tools"
docker save -o "$output_folder"/bioinfo_tools.tar bioinfo_tools

docker build -t custom_sra_tools "$repository_path/dockerfiles/custom_sra_tools"
docker save -o "$output_folder"/custom_sra_tools.tar custom_sra_tools
