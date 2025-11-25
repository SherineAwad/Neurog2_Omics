# Neurog2 Multiomics Using Seurat and Signac

```text
 _   _                            ____  
| \ | | ___ _   _ _ __ ___   __ _|___ \ 
|  \| |/ _ \ | | | '__/ _ \ / _` | __) |
| |\  |  __/ |_| | | | (_) | (_| |/ __/ 
|_| \_|\___|\__,_|_|  \___/ \__, |_____|
                            |___/       
```

Normally, adult mouse retinas cannot generate new neurons, so vision loss is permanent. We are testing whether overexpressing Neurog2 in Müller glia can reprogram them into functional retinal neurons. These newly formed neurons behave like normal retinal cells.

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


## Remove cluster 11, 19, and 26 and reCluster 


![reCluster by sample](figures/cNeurog2_Recluster_BySample.png?v=7)

![reCluster by Clusters](figures/cNeurog2_Recluster_Clusters.png?v=6)


### UMI per cluster for investigating poor cells


![](figures/cNeurog2_UMI_Violin_HarmonyWNN.png?v=3)

## Marker Genes feature  plots


<img src="figures/cNeurog2_Abca8a_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Acta2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Apoe_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Aqp4_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Arr3_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Ascl1_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Atoh7_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Bsn_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Cabp5_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Calb1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Calb2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Ccr2_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Chat_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Csf1r_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Elavl4_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Emx1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Gad1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Gfap_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Glul_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Hes1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Hes5_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Insm1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Isl1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Lhx1_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Lhx2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Lhx4_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Malat1_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_mScarlet3_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_mt-Atp6_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Neurog2-9SA_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Neurog2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Notch1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Nrl_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Olig2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Otx2_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Pax2_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Pax6_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Prdm1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Prdx6_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Rbfox3_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Rho_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Rlbp1_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Sebox_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Slc17a7_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Slc1a3_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Slc6a9_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Sox11_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Sox9_FeaturePlot.png?v=29" width="200">

<img src="figures/cNeurog2_Tfap2a_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Tie1_FeaturePlot.png?v=29" width="200"> <img src="figures/cNeurog2_Vim_FeaturePlot.png?v=29" width="200">



## Annotation: 


![Annotations](figures/cNeurog2_annotated_umap.png?v=9)

## Cell ratio  

![Cell ratio](figures/cNeurog2_celltype_ratio_by_sample.png?v=9)

## Differential Gene Expression

#### Starting with: 

| Metric          | Count   |
|-----------------|---------|
| Number of cells | 17,278  |
| Number of genes | 32,287  |



```r
celltype_markers <- FindAllMarkers(
  myObject,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

``` 

# 📝📝📝📝 Download DGE 

#### Ascl1 exist chat not 

[Download DGE csv](https://docs.google.com/spreadsheets/d/14xipS-nTOasxoGZ4Ljx0jYxBgCW2svklyGPC1oQ7KPQ/edit?usp=sharing)

![](figures/Top10_Markers_heatmap.png?v=5)


### Heatmap using `hPlot.R` script

![Top 50 genes heatmap](figures/Top50_pval.rna.Avg_SCT_heatmap_Zscore.png?v=2)

# Differential Peaks

```
myObject.atac.markers <- FindAllMarkers(
  myObject,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
``` 

#### Download unfiltered annotated Diffpeaks 

[All Annotated DiffPeaks](https://docs.google.com/spreadsheets/d/1SM759_168C0F-RcrduTf4zr1dbh4U2Ta_WPNwRQ0YlE/edit?usp=sharing)


![](figures/TopPeaksNearbyheatmap.png?v=1)


#  Subset MG/MGPC 

## Differential gene expression in the MG/MGPC subset 

```r
dge_results <- FindMarkers(
  subset_cells,
  ident.1 = "TH2",
  ident.2 = "TH1",
  assay = "SCT",
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  )

``` 

![Top markers heatmap](figures/cNeurog2_TH2_vs_TH1_heatmap_top30.png?v=4) 


## Download DGE for MG/MGPC subset 

[MG/MGPC DGE](https://docs.google.com/spreadsheets/d/17JkCb6IviAh6rUSZlV8lAQTN0331oyK8BmUCkgzcvnA/edit?usp=sharing)



## Differential Peaks in subset MG/MGPC 


## Download annotated MG/MGPC subset differential peaks

[MG/MGPC differential peaks](https://docs.google.com/spreadsheets/d/17LYiDVUW5CcdYlOw8LkqRCJfDMJGFKz0Wao67wdDOqE/edit?usp=sharing)


Neurog2 one close peak which has the following significance 
 
| Peak                     | p_val        | avg_log2FC | pct.1 | pct.2 | p_val_adj   |
| ------------------------ | ------------ | ---------- | ----- | ----- | ----------- |
| chr3-127633860-127634766 | 8.70 × 10⁻¹⁵ | -0.4201    | 0.135 | 0.193 | 1.85 × 10⁻⁹ |


However, its ranked 1+k in order !!



## MG/MGPC Differential Peaks heatmap 

![MG/MGPC diffPeaks](figures/cNeurog2_MG_MGPC_diffPeaksheatmap.png?v=3)
 

## Motif enrichment on MG/MGPC subset 

```
**Running motif enrichment on UPREGULATED peaks...**  
**Valid UPREGULATED peaks with GC content:** 0  
❌ **Error:** No valid UPREGULATED peaks with GC content data


The GC content error show up even with Upregulated diff peaks, to fix this we do, skip GC check: 

Running motif enrichment on UPREGULATED peaks...
Warning: No GC content data available. Running motif enrichment without GC correction.
Selecting background regions to match input sequence characteristics
Matching GC.percent distribution
Testing motif enrichment in 2279 regions
```


![UP MG/MGPC Motifs](figures/MG_MGPC_UpMotifEnrich.png?=v1)

# 📝📝📝📝 Download MG/MGPC Upregulated Motifs 

[UP MG/MGPC Motifs](https://docs.google.com/spreadsheets/d/1Oxw1cfRA5kpBqGrrQw34i3wvBQvyKzwDmzUj1YuCY3U/edit?usp=sharing)
 
### 🚨 🚨 🚨 Warning:  Neurog2 is not still found 

# Pathways of genes near upregulated peaks in TH2 

![GO](figures/TH2_pathways_dotplot.png)

# More plots ammendment 

```
Rscript hPlot.R myobject.rds topgenes.txt 
```

![]()

# Coverage plot for some genes in MG/MGPC Subset Control/Overexpressed

###### Neurog2 
![Neuorg2](figures/Neurog2_MG_MGPC_coverage.png)

##### Rho 
![Rho](figures/Rho_MG_MGPC_coverage.png)


##### Dlx1 
![Dlx1](figures/Dlx1_MG_MGPC_coverage.png)

#### Opn1sw
![Opn1sw](figures/Opn1sw_MG_MGPC_coverage.png)


# Coverage plot some genes in all cell types

##### Neurog2 
![Neurog2_ALL](figures/Neurog2_coverage.png)


