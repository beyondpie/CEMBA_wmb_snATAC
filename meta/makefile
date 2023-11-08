SHELL = /bin/bash
.PHONY = clean install

install: CEMBA_unanalyzed_samples.txt

## CEMBA_analyzed.tab from
## /projects/ps-renlab/yangli/projects/CEMBA/01.cluster/mba.snap.list
## which is the samples analyzed by Yang.

## old.sample.ls from
## /projects/ps-renlab/yangli/projects/CEMBA/01.cluster/old.sample.ls
## which is the CEMBA v1 samples.

CEMBA_analyzed_samples.txt: CEMBA_analyzed.tab
	awk -F "\t" '{print $$1}' $< | sort > $@

CEMBA_tracking_samples.txt: CEMBA_tracking.csv
	-dos2unix -n $< CEMBA_tracking_unix.csv
	awk -F ',' '{print $$2}' CEMBA_tracking_unix.csv | sort > $@

old_samples.txt: old.sample.list
	sort $< > $@

CEMBA_unanalyzed_samples.txt: CEMBA_tracking_samples.txt CEMBA_analyzed_samples.txt old_samples.txt
	comm -23 CEMBA_tracking_samples.txt CEMBA_analyzed_samples.txt > CEMBA_unanalyzed_samples_from_Yang.txt && \
	comm -23 CEMBA_tracking_samples.txt old_samples.txt > CEMBA_not_old_samples.txt && \
  comm -12 CEMBA_unanalyzed_samples_from_Yang.txt CEMBA_not_old_samples.txt > $@
	- rm CEMBA_unanalyzed_samples_from_Yang.txt
	- rm CEMBA_not_old_samples.txt

gtf := modified_gencode.vM23.primary_assembly.annotation.gtf
bedfile := mouse.modified.gencode.vM23.bed
geneUp2kfile := mouse.modified.gencode.vM23.gene.up2k.bed

${bedfile}:${gtf}
	awk 'BEGIN{FS=OFS="\t"}($$3=="gene"){split($$9,a,"\""); print $$1,$$4-1,$$5,a[8],$$6,$$7}' $< \
    | sort -k1,1 -k2,2n > $@

${geneUp2kfile}:${bedfile}
	awk 'BEGIN{FS=OFS="\t"}{if($$6=="+" && $$2-2000>0){print $$1,$$2-2000,$$3,$$4,$$5,$$6}else if($$6=="+" && $$2-2000<0){print $$1,0,$$3,$$4,$$5,$$6}else if($$6=="-"){print $$1,$$2,$$3+2000,$$4,$$5,$$6}}' $< > $@

.PHONY: geneUp2k
geneUp2k:${geneUp2kfile}

bedfile2 := mouse.modified.gencode.vM23.ENSMUSG.bed
geneUp2kfile2 := mouse.modified.gencode.vM23.ENSMUSG.up2k.bed
${bedfile2} : ${gtf}
	awk 'BEGIN{FS=OFS="\t"}($$3=="gene"){split($$9,a,"\""); print $$1,$$4-1,$$5,a[2],$$6,$$7}' $< \
    | sort -k1,1 -k2,2n > $@
${geneUp2kfile2}:${bedfile2}
	awk 'BEGIN{FS=OFS="\t"}{if($$6=="+" && $$2-2000>0){print $$1,$$2-2000,$$3,$$4,$$5,$$6}else if($$6=="+" && $$2-2000<0){print $$1,0,$$3,$$4,$$5,$$6}else if($$6=="-"){print $$1,$$2,$$3+2000,$$4,$$5,$$6}}' $< > $@
.PHONY: geneUp2k2
geneUp2k2:${geneUp2kfile2}

ensymbol2symbol_vM23.tsv: ${gtf}
	awk 'BEGIN{FS=OFS="\t"}{if($$3=="gene"){split($$9, a,"\""); print a[2],a[8]}}' $< > $@

# used for igv to load gtf
# https://www.biostars.org/p/289819/
mm10.modified.vM23.sorted.gtf.gz: modified.genecode.vM23.primary_assembly.annotation.gtf
	(grep ^"#" $<; grep -v ^"#" $< | sort -k1,1 -k4,4n) | bgzip  > $@
	tabix -p gff $@

## update vm23
# sed 's/;$//' modified_gencode.vM23.primary_assembly.annotation.gtf > gencode.vM23.primary.rm.semicolon.gtf
# awk -F '\t' '{if( ($3 ~ /gene/) || (NR < 6)) {print $0}}' gencode.vM23.primary.rm.semicolon.gtf  > gencode.vM23.gene.annot2.rm.semicolon.gtf 
# awk -F '\t' '{if( $0 ~ /intro/) {print $0";"} else {print $0}}' modified_gencode.vM23.primary_assembly.annotation.gtf > gencode.vM23.primary.annot2.gtf 

clean:
	-rm old_samples.txt
	-rm CEMBA_analyzed_samples.txt
	-rm CEMBA_tracking_samples.txt
	-rm CEMBA_tracking_unix.csv

