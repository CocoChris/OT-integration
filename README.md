# OT-integreation
Integration of multiple single-cell datasets based on optimal transport

## pbmc
- run cca_pbmc.R to determine the optimal number of CCA vectors and get CCA results
- run select_optimal_cluster_number.R to select the optimal clustering number of 10X dataset
- run calcBIC.R to determine the optimal clustering number of the mean partition, the number is also treated as the final clustering number in the Seurat and LIGER workflow
- run seurat3_workflow_pbmc.R to get tsne and clustering results based on Seurat
- run liger_workflow_pbmc.R to get tsne and clustering results based on LIGER
- run integration_pbmc.m to integrate the datasets based on OT and get the gamma matrix
- run makeRiverPlot_pbmc.R to draw the riverplot on the gamma matrix

## human pancreas
- run seurat_workflow_human_pancreas.R to get cca, tsne and clustering results based on Seurat
- run liger_workflow_human_pancreas.R to get tsne and clustering results based on LIGER
- run integration_human_pancreas.m to integrate the datasets based on OT and get the gamma matrix
- run makeRiverPlot_human_pancreas.R to draw the riverplot on the gamma matrix
