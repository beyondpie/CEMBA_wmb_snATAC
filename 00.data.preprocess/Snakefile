##############################################
#global settings
shell.prefix("set -o pipefail; ")
shell.prefix("set -e; ")
shell.prefix("set -u; ")
localrules: all
configfile: "config.json"

# link
from snakemake.exceptions import MissingInputException

#include_prefix = "/projects/ps-renlab/yangli/scripts/snakemake/"

#include:
#    include_prefix + "rules/"

dissectionItems = {
    datasets: dissections for datasets, dissections in config["datasets"].items()
}

#print(dissectionItems.keys())
#print(dissectionItems.values())

rule all:
    input:
#        expand("{dissection}/{dataset}/{dataset}.snap.cluster.RData", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
#        expand("{dissection}/{dataset}/{dataset}.snap.plotGene.pdf", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
#        expand("{dissection}/{dataset}/{dataset}.snap.bed.gz", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
#        expand("{dissection}/{dataset}/{dataset}.snap.add_bmat.ok", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
#        expand("{dissection}/{dataset}/{dataset}.snap.add_gmat.ok", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
        expand("{dissection}/{dataset}/processed/{dataset}.bam", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values()),
        expand("{dissection}/{dataset}/processed/{dataset}.bedpe.gz", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values())
#        expand("{dissection}/{dataset}/snap2cb/{dataset}.snap.exprMatrix.tsv.gz", zip, dataset=dissectionItems.keys(), dissection=dissectionItems.values())


rule snap_align: 
    input:
        fq1 = "{dissection}/{dataset}/raw/{dataset}.demultiplexed.R1.fastq.gz",
        fq2 = "{dissection}/{dataset}/raw/{dataset}.demultiplexed.R2.fastq.gz"
    output:
        bam = "{dissection}/{dataset}/processed/{dataset}.bam"
    params:
        ref = config["references"]["REF_FA"],
        aligner = config["snap_align"]["aligner"],
        path_to_aligner = config["snap_align"]["path_to_aligner"],
        read_fastq_command = config["snap_align"]["read_fastq_command"],
        min_cov = config["snap_align"]["min_cov"],
        num_threads = config["snap_align"]["num_threads"],
        if_sort = config["snap_align"]["if_sort"],
        tmp_folder = config["snap_align"]["tmp_folder"],
        overwrite = config["snap_align"]["overwrite"],
        jobname = "{dataset}.align"
    benchmark:
        "benchmarks/{dissection}/{dataset}.align.benchmark"
    log:
        "log/{dissection}/{dataset}.align.log"
    shell:
        """
        /home/yangli1/apps/anaconda2/bin/snaptools align-paired-end  \
        --input-reference={params.ref}  \
        --input-fastq1={input.fq1} \
        --input-fastq2={input.fq2}  \
        --output-bam={output.bam}  \
        --aligner={params.aligner}  \
        --path-to-aligner={params.path_to_aligner}  \
        --read-fastq-command={params.read_fastq_command}  \
        --min-cov={params.min_cov}  \
        --num-threads={params.num_threads}  \
        --if-sort={params.if_sort}  \
        --tmp-folder={params.tmp_folder}  \
        --overwrite={params.overwrite} &> {log}
        """

rule snap_pre:
    input:
        bam="{dissection}/{dataset}/processed/{dataset}.bam"
    output:
        snap="{dissection}/{dataset}/processed/{dataset}.snap"
    params:
        chrom_sizes=config["references"]["REF_CHROM_SIZE"],
        jobname = "{dataset}.snap_pre"
    benchmark:
        "benchmarks/{dissection}/{dataset}.snap_pre.benchmark"
    log:
        "log/{dissection}/{dataset}.snap_pre.log"
    shell:
        """
        /home/yangli1/apps/anaconda2/bin/snaptools snap-pre  \
        --input-file={input.bam}  \
        --output-snap={output.snap} \
        --genome-name=mm10  \
        --genome-size={params.chrom_sizes}  \
        --min-mapq=30  \
        --min-flen=0  \
        --max-flen=1000  \
        --keep-chrm=TRUE  \
        --keep-single=TRUE  \
        --keep-secondary=False  \
        --overwrite=True  \
        --min-cov=100  \
        --verbose=True &> {log}
        """

rule snap_add_bmat:
    input:
        snap="{dissection}/{dataset}/processed/{dataset}.snap"
    output:
        "{dissection}/{dataset}/processed/{dataset}.snap.add_bmat.ok"
    params:
        jobname = "{dataset}.snap_add_bmat"
    benchmark:
        "benchmarks/{dissection}/{dataset}.snap_add_bmat.benchmark"
    log:
        "log/{dissection}/{dataset}.snap_add_bmat.log"
    shell:
        """
        /home/yangli1/apps/anaconda2/bin/snaptools snap-add-bmat \
        --snap-file={input.snap} \
        --bin-size-list 1000 5000 10000 \
        --verbose=True &> {log}
        echo "Done" > {output}
        """

rule snap_add_gmat:
    input:
        snap="{dissection}/{dataset}/processed/{dataset}.snap",
        bmatok="{dissection}/{dataset}/processed/{dataset}.snap.add_bmat.ok"
    output:
        "{dissection}/{dataset}/processed/{dataset}.snap.add_gmat.ok"
    params:
        anno_bed=config["references"]["REF_ANNO_BED"],
        jobname = "{dataset}.snap_add_gmat"
    benchmark:
        "benchmarks/{dissection}/{dataset}.snap_add_gmat.benchmark"
    log:
        "log/{dissection}/{dataset}.snap_add_gmat.log"
    shell:
        """
        /home/yangli1/apps/anaconda2/bin/snaptools snap-add-gmat \
        --snap-file={input.snap} \
        --gene-file={params.anno_bed} \
        --verbose=True &> {log}
        echo "Done" > {output}
        """

rule pre_sta:
    input:
#        bmatok="{dissection}/{dataset}/{dataset}.snap.add_bmat.ok",
#        gmatok="{dissection}/{dataset}/{dataset}.snap.add_gmat.ok",
        snap="{dissection}/{dataset}/processed/{dataset}.snap"
    output:
        RData="{dissection}/{dataset}/processed/{dataset}.snap.pre.RData",
        summary="{dissection}/{dataset}/processed/{dataset}.snap.qc.summary"
    params:
        snapATACpre=config["scripts"]["snapATACpre"],
        blacklist=config["references"]["REF_BLACKLIST"],
        cpus=config["pre_sta"]["cpu"],
        jobname = "{dataset}.pre_sta"
    benchmark:
        "benchmarks/{dissection}/{dataset}.pre_sta.benchmark"
    log:
        "log/{dissection}/{dataset}.pre_sta.log"
    shell:
        """
        Rscript {params.snapATACpre} -i {input.snap} --fragment_num 1000 --umi_num 500 --mito_ratio 1 --dup_ratio 1 --umap_ratio 0 --pair_ratio 0 --bin_size 5000 --black_list {params.blacklist} --cpu {params.cpus} -o {input.snap} 1> {log} 2> {output.summary}
        """

rule cluster:
    input:
        RData="{dissection}/{dataset}/processed/{dataset}.snap.pre.RData"
    output:
        cluster="{dissection}/{dataset}/processed/{dataset}.snap.cluster.RData"
    params:
        prefix="{dissection}/{dataset}/processed/{dataset}.snap",
        snapATACcluster=config["scripts"]["snapATACcluster"],
        jobname = "{dataset}.cluster"
    benchmark:
        "benchmarks/{dissection}/{dataset}.cluster.benchmark"
    log:
        "log/{dissection}/{dataset}.cluster.log"
    shell:
        """
        Rscript {params.snapATACcluster} -i {input.RData} -o {params.prefix} &> {log}
        """

rule plotGene:
    input:
        RData="{dissection}/{dataset}/processed/{dataset}.snap.cluster.RData"
    output:
        plotGene="{dissection}/{dataset}/processed/{dataset}.snap.plotGene.pdf"
    params:
        prefix="{dissection}/{dataset}/processed/{dataset}.snap",
        snapATACplotGene=config["scripts"]["snapATACplotGene"],
        markerGenes=config["plotGene"]["markerGenes"],
        jobname = "{dataset}.plotGene"
    benchmark:
        "benchmarks/{dissection}/{dataset}.plotGene.benchmark"
    log:
        "log/{dissection}/{dataset}.plotGene.log"
    shell:
        """
        Rscript {params.snapATACplotGene} -i {input.RData} -m {params.markerGenes} -o {params.prefix} &> {log}
        """

rule dump_frag:
    input:
        bmatok="{dissection}/{dataset}/processed/{dataset}.snap.add_bmat.ok",
        gmatok="{dissection}/{dataset}/processed/{dataset}.snap.add_gmat.ok",
        snap="{dissection}/{dataset}/processed/{dataset}.snap"
    output:
        bed="{dissection}/{dataset}/processed/{dataset}.snap.bed.gz"
    params:
        jobname = "{dataset}.dump_frag"
    benchmark:
        "benchmarks/{dissection}/{dataset}.dump_frag.benchmark"
    log:
        "log/{dissection}/{dataset}.dump_frag.log"
    shell:
        """
        /home/yangli1/apps/anaconda2/bin/snaptools dump-fragment --snap-file {input.snap} --output-file {output.bed} &> {log}
        """

rule tsse2depth:
    input:
        bam = ancient("{dissection}/{dataset}/processed/{dataset}.bam")
    output:
        ok = "{dissection}/{dataset}/processed/{dataset}.tsse2depth.ok",
        bam = "{dissection}/{dataset}/processed/filtered.bam"
    params:
        outdir = "{dissection}/{dataset}/processed/",
        jobname = "{dataset}.tsse"
    benchmark:
        "benchmarks/{dissection}/{dataset}.tsse2depth"
    log:
        "log/{dissection}/{dataset}.tsse2depth.log"
    shell:
        """
        /projects/ps-renlab/yangli/scripts/snATACutils/bin/calTSSePE {input.bam} {params.outdir} --pair -g /projects/ps-renlab/yangli/genome/mm10/gencode.vM16.annotation.gtf
        echo 'Done' > {output.ok}
        """


rule bam2bedpe:
    input:
        bam = "{dissection}/{dataset}/processed/filtered.bam"
    output:
        bedpe = "{dissection}/{dataset}/processed/{dataset}.bedpe.gz"
    params:
        jobname = "{dataset}.bam2bedpe",
        cleanfix = "{dissection}/{dataset}/processed/filtered.nsrt.cleanfix.bam"
    benchmark:
        "benchmarks/{dissection}/{dataset}.bam2bedpe.benchmark"
    log:
        "log/{dissection}/{dataset}.bam2bedpe.log"
    shell:
        """
        samtools sort -n {input.bam} | samtools view -b -F 2048 -f 3 - | samtools fixmate -r -m - {params.cleanfix}
        bedtools bamtobed -bedpe -mate1 -i {params.cleanfix} | gzip > {output.bedpe}
        """

rule convertBed:
    input:
        bed="{dissection}/{dataset}/processed/{dataset}.snap.bed.gz"
    output:
        bam="{dissection}/{dataset}/processed/{dataset}.snap.srt.bam",
        bw="{dissection}/{dataset}/processed/{dataset}.snap.bw"
    params:
        chrom_sizes=config["references"]["REF_CHROM_SIZE"],
        jobname = "{dataset}.convertBed"
    benchmark:
        "benchmarks/{dissection}/{dataset}.convertBed.benchmark"
    log:
        "log/{dissection}/{dataset}.convertBed.log"
    shell:
        """
        zcat {input.bed} | bedToBam -i - -g {params.chrom_sizes} | samtools sort - > {output.bam}
        samtools index -b {output.bam}
        bamCoverage --bam {output.bam} -o {output.bw} \
        --binSize 10
        --normalizeUsing RPGC
        --effectiveGenomeSize 2150570000
        --ignoreForNormalization chrX chrM
        --extendReads
        """

rule snap2cb:
    input:
        Rdata=ancient("{dissection}/{dataset}/processed/{dataset}.snap.gmat.RData")
    output:
        exprMtx="{dissection}/{dataset}/processed/snap2cb/{dataset}.snap.exprMatrix.tsv.gz"
    params:
        dat="{dataset}",
        prefix="{dissection}/{dataset}/processed/snap2cb/{dataset}.snap",
        geneF="{dissection}/{dataset}/processed/snap2cb/{dataset}.snap.gene.tsv",
        cellF="{dissection}/{dataset}/processed/snap2cb/{dataset}.snap.cell.tsv",
        mtxF="{dissection}/{dataset}/processed/snap2cb/{dataset}.snap.gmat.mtx",
        metaF="CEMBA.data.201906.summary",
        summaryHtml="{dissection}/{dataset}/processed/snap2cb/summary.html",
        methodsHtml="{dissection}/{dataset}/processed/snap2cb/methods.html",
        downloadsHtml="{dissection}/{dataset}/processed/snap2cb/downloads.html",
        snapATACsnap2cb=config["scripts"]["snapATACsnap2cb"],
        jobname = "{dataset}.snap2cb"
    benchmark:
        "benchmarks/{dissection}/{dataset}.snap2cb.benchmark"
    log:
        "log/{dissection}/{dataset}.snap2cb.log"
    shell:
        """
        Rscript {params.snapATACsnap2cb} -i {input.Rdata} -m {params.metaF} -o {params.prefix}
        cbTool mtx2tsv {params.mtxF} {params.geneF} {params.cellF} {output.exprMtx}
        echo -e "CEMBA project\\n\\nThis result is based on: " {params.dat} > {params.summaryHtml}
        echo 'snATAC-seq: Combinatorial Barcoding' > {params.methodsHtml}
        echo 'track files: http://renlab.sdsc.edu/yangli/cbHub/' > {params.downloadsHtml}
        echo 'Done'
        """



