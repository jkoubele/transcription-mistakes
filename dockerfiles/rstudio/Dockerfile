FROM rocker/rstudio:4.4.1
RUN Rscript -e 'install.packages("tidyverse")'
RUN apt-get update
RUN apt-get install libglpk-dev -y
RUN Rscript -e 'install.packages("Seurat")'
RUN Rscript -e 'install.packages("jsonlite")'