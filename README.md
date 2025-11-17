# Neurog2 Multiomics Using Seurat and Signac

Overexpression of stablized version of Neurog2 


## Before filtering 

![TH1 before filtering](figures/TH1_QC_vlnplot.png?v=4)


![TH2 before filtering](figures/TH2_QC_vlnplot.png?v=4)



## After filtering 

![TH1 after filtering](figures/TH1_QC_vlnplot_after_filtering.png?v=4)


![TH2 after filtering](figures/TH2_QC_vlnplot_after_filtering.png?v=4)


## Umap 

![Clusters Vln1](figures/cNeurog2_Cluster_VlnPlot1.png?v=4)

![clusters Vln2](figures/cNeurog2_Cluster_VlnPlot2.png?v=4)

![Clusters](figures/cNeurog2_Clusters.png?v=4)

![WNN by Sample](figures/cNeurog2_WNN_by_sample.png?v=3)


## Samples are too seperated, we do harmony for batch effect 

![Clusters harmony Vln1](figures/cNeurog2_Cluster_VlnPlot2_harmony.png?v=3)

![Clusters harmony Vln2](figures/cNeurog2_Cluster_VlnPlot1_harmony.png?v=3)

![Clusters harmony clusters](figures/cNeurog2_Clusters_harmony.png?v=3)

![Clusters harmony by sample](figures/cNeurog2_WNN_by_sample_harmony.png?v=3)


## Remove cluster 11, 19, and 20 and reCluster 


![reCluster by sample](figures/cNeurog2_Recluster_BySample.png?v=4)

![reCluster by Clusters](figures/cNeurog2_Recluster_Clusters.png?v=4)



### UMI per cluster for investigating poor cells


![](figures/cNeurog2_UMI_Violin_HarmonyWNN.png?v=1)

## Marker Genes feature  plots


<img src="figures/cNeurog2_Abca8a_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Acta2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Apoe_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Aqp4_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Arr3_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Ascl1_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Atoh7_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Bsn_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Cabp5_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Calb1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Calb2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Ccr2_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Chat_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Csf1r_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Elavl4_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Emx1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Gad1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Gfap_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Glul_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Hes1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Hes5_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Insm1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Isl1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Lhx1_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Lhx2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Lhx4_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Malat1_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_mScarlet3_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_mt-Atp6_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Neurog2-9SA_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Neurog2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Notch1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Nrl_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Olig2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Otx2_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Pax2_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Pax6_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Prdm1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Prdx6_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Rbfox3_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Rho_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Rlbp1_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Sebox_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Slc17a7_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Slc1a3_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Slc6a9_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Sox11_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Sox9_FeaturePlot.png?v=25" width="200">

<img src="figures/cNeurog2_Tfap2a_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Tie1_FeaturePlot.png?v=25" width="200"> <img src="figures/cNeurog2_Vim_FeaturePlot.png?v=25" width="200">



## Annotation: 


![Annotations](figures/cNeurog2_annotated_umap.png?v=4)

## Cell ratio  

![Cell ratio](figures/cNeurog2_celltype_ratio_by_sample.png?v=4)

## Differential Gene Expression

![](figures/cNeurog2_celltype_markers_heatmap.png?v=3)

![](figures/cNeurog2_celltype_top_markers_featureplot.png?v=3)

![](figures/cNeurog2_celltype_markers_dotplot.png?v=3)

##### Rod 
![](figures/cNeurog2_Rod_markers_violin.png?v=3)

#### BC

![](figures/cNeurog2_BC_markers_violin.png?v=3)

#### MGPC
![](figures/cNeurog2_MGPC_markers_violin.png?v=3)

#### MG
![](figures/cNeurog2_MG_markers_violin.png?v=3)


## Top Differential Peaks after annotation (nearby Genes)

#### Starting with: 

| Metric          | Count   |
|-----------------|---------|
| Number of cells | 17,360  |
| Number of genes | 213,119 |


![](figures/TopPeaks.png?v=2)

![](figures/peak_category_by_cluster.png?v=1)

![](figures/peak_categories.png?v=1)
