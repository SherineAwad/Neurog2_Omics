configfile: "config.yaml"

SAMPLES = config["samples"]
MERGE = config["merge_name"]
P = config["filter_params"]

rule all:
    input:
        config["diffpeaks_annotated"],
        config["peak_plot"]

# --------------------------
# Preprocess each sample
# --------------------------
rule preprocess:
    output:
        touch("{sample}.preprocessed")
    shell:
        "Rscript preprocessS.R {wildcards.sample} && touch {output}"

# --------------------------
# Filter each sample
# --------------------------
rule filter:
    input:
        "{sample}.preprocessed"
    output:
        touch("{sample}.filtered")
    params:
        nATAC_max=P["nATAC_max"],
        nRNA_max=P["nRNA_max"],
        nATAC_min=P["nATAC_min"],
        nRNA_min=P["nRNA_min"],
        nFeature_RNA_min=P["nFeature_RNA_min"],
        nucleosome_signal_max=P["nucleosome_signal_max"],
        percent_mt_max=P["percent_mt_max"],
        tss=P["tss"]
    shell:
        ("Rscript filterS.R {wildcards.sample} "
         "{params[nATAC_max]} {params[nRNA_max]} "
         "{params[nATAC_min]} {params[nRNA_min]} "
         "{params[nFeature_RNA_min]} {params[nucleosome_signal_max]} "
         "{params[percent_mt_max]} {params[tss]} && touch {output}")

# --------------------------
# Merge samples
# --------------------------
rule merge:
    input:
        expand("{sample}.filtered", sample=SAMPLES)
    output:
        touch(f"{MERGE}.merged")
    shell:
        "Rscript mergeS.R {MERGE} {SAMPLES[0]} {SAMPLES[1]} && touch {output}"

# --------------------------
# Analyze merged
# --------------------------
rule analyse:
    input:
        f"{MERGE}.merged"
    output:
        touch(f"{MERGE}.analysed")
    shell:
        "Rscript analyseS.R {MERGE} && touch {output}"

# --------------------------
# Recluster
# --------------------------
rule recluster:
    input:
        f"{MERGE}.analysed"
    output:
        touch(f"{MERGE}.reclustered")
    shell:
        "Rscript reCluster.R {MERGE} && touch {output}"

# --------------------------
# Plot markers
# --------------------------
rule plot_markers:
    input:
        f"{MERGE}.reclustered"
    output:
        touch(f"{MERGE}.markers_plotted")
    params:
        markers=config["markers_file"]
    shell:
        "Rscript plotMarkersS.R {MERGE} {params.markers} && touch {output}"

# --------------------------
# Annotate
# --------------------------
rule annotate:
    input:
        f"{MERGE}.reclustered"
    output:
        touch(f"{MERGE}.annotated")
    params:
        annotations=config["annotations_file"]
    shell:
        "Rscript annotate.R {MERGE} {params.annotations} && touch {output}"

# --------------------------
# Differential gene expression
# --------------------------
rule dge:
    input:
        f"{MERGE}.reclustered"
    output:
        touch(f"{MERGE}.dge")
    shell:
        "Rscript dge.R {MERGE} && touch {output}"

# --------------------------
# Differential peaks
# --------------------------
rule diff_peaks:
    input:
        f"{MERGE}.reclustered"
    output:
        config["diffpeaks_file"]
    shell:
        "Rscript diffPeaks.R {MERGE}"

# --------------------------
# Annotate peaks
# --------------------------
rule annotate_peaks:
    input:
        config["diffpeaks_file"],
        config["gtf_file"]
    output:
        config["diffpeaks_annotated"]
    shell:
        "Rscript annotatePeaks.R --file {input[0]} --gtf {input[1]} --out {output}"

# --------------------------
# Plot peaks
# --------------------------
rule plot_peaks:
    input:
        config["diffpeaks_annotated"]
    output:
        config["peak_plot"]
    shell:
        "Rscript plotPeaks.R -f {input} -n 5 -o {output}"

