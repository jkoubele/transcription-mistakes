library(Seurat)
library(tidyverse)
library(jsonlite)

gse_list <- c('GSE147319', 'GSE158859', 'GSE190848', 'GSE230816')

mmus_cell_cycle_genes <- fromJSON("/data/mmus_cell_cycle_genes.json")
for (gse in gse_list){
  print(gse)
  gsm_list <- list.dirs(path = paste0("/data/aligned/",gse), full.names = FALSE, recursive = FALSE)
  for (gsm in gsm_list){
    star_solo_path <- paste0("/data/aligned/", gse, "/", gsm, "/Solo.out/GeneFull/filtered")
    print(star_solo_path)
    seurat_data <- ReadSTARsolo(star_solo_path)
    seurat_object <- CreateSeuratObject(counts = seurat_data)

    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
    VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    seurat_object <- NormalizeData(seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
    seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
    seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))

    seurat_object <- CellCycleScoring(seurat_object,
                                      s.features = mmus_cell_cycle_genes$mmus_s,
                                      g2m.features = mmus_cell_cycle_genes$mmus_g2m,
                                      set.ident = TRUE)

    output_folder <- paste0("/data/cell_cycle_info/", gse, "/", gsm)
    dir.create(output_folder, recursive = TRUE)
    write_delim(seurat_object@meta.data, paste0(output_folder, "/metadata.tsv"), delim = '\t')

  }
}
