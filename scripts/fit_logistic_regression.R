library(tidyverse)

X <- read.delim('/cellfile/datapublic/jkoubele/sc_data/design_matrices/GSE190848/design_matrix.tsv')
print(paste("Average error rate:",weighted.mean(X$y, X$count) * 1e6))

model <- glm(y~base+base_next+base_prev+phase+gsm, weights=count, data=X, family = "binomial")
# model <- glm(y~1, weights=count, data=X, family = "binomial")
summary(model)