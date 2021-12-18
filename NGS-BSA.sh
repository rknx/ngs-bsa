#!/bin/sh

echo -e '\nNGS-BSA resistance gene mapping by Anuj Sharma.\n'

########## Preparing input parameters ##########

params()
{
   echo ''
   echo -e 'Usage:\n'$0' [-h <help>] -g gene -s species [-o output folder] -t resistance type -b {bulk codes} {[-l lane names]} [-r reference folder] -v reference line [-f reference genome link] [-a annotation link] [-e snpEff config folder] [-k known snps link] [-c cut-off threshold] [-w smoothing span] [-i iteration] [-x steps to skip] [-d steps to redo] [-n <nolabel>] [-p path to script folder]'
   echo -e '\nwhere:'
   echo -e '\t-g\tShort name of gene under analysis (alphanumeric)'
   echo -e '\t-s\tSpecies name, (alphanumeric and "_")'
   echo -e '\t-o\tOutput folder (alphanumeric and "_")'
   echo -e '\t-t\tResistance type ("d: dominant", "r: recessive" or "c: codominant")'
   echo -e '\t-b\tList of names or codes for bulks or samples (alphanumeric and "_").\n\t\tInput resitant bulk followed by susceptible bulk.\n\t\tOnly first two bulks will be used for plotting.'
   echo -e '\t-l\tLane names if done on multiple lanes (alphanumeric and "_")'
   echo -e '\t-r\tReference directory path (Created if not present)'
   echo -e '\t-v\tReference variety (alphanumeric and "_")'
   echo -e '\t-f\tReference genome fasta gzip download link'
   echo -e '\t-a\tAnnotation file gzip link for reference genome'
   echo -e '\t-e\tsnpEff config file path'
   echo -e '\t-k\tKnown SNPs vcf gzip download link'
   echo -e '\t-c\tCut-off threshold for manhattan plot (numeric, 0-1, default = 0.3)'
   echo -e '\t-w\tSpanning width for smoothing in Manhattan plot (numeric, 0-1, default = 0.2)'
   echo -e '\t-i\tRun iteration (numeric, default=NULL)'
   echo -e '\t-x\tManually skip some steps. Run -x help for list of steps'
   echo -e '\t-d\tForce redo some steps. Run -d help for list of steps'
   echo -e '\t-n\tRemove gene label in plots (numeric, default=NULL)'
   echo -e '\t-p\tPath to script folder'
   echo ''
   echo -e 'Instructions:\n\tInstall the following programs: Java, Samtools, BWA, SnpEff, GATK and R.\n\tInstall following packages in R: ggplot2, ggrepel, reshape2.\n\tPlace fastq files (fastq.gz, fastq, fq.gz or fq) from both pools in root folder.\n\tRun script with appropriate parameters.\n\tIf using SLURM, fill in the variables in batch file and run this script from batch file.'
   echo ''
}

export output="./output"
refdir="./refs"
snpeffdir="./snpEff"
export lanes=()
thres=0.3
iteration=1
label="y"
span=0.2

while getopts hg:s:o:t:b:l:r:v:f:a:e:k:c:w:i:x:d:np: args #jmquyz
do
case "${args}" in
h) params;;
g) export gene="$OPTARG";;
s) species="$OPTARG";;
o) export output="$OPTARG";;
t) rtype="$OPTARG";;
b) bulks+=("$OPTARG");;
l) lanes+=("$OPTARG");;
r) refdir="$OPTARG";;
v) export refline="$OPTARG";;
f) refdl="$OPTARG";;
a) annotdl="$OPTARG";;
e) snpeffdir="$OPTARG";;
k) snpdl="$OPTARG";;
c) thres="$OPTARG";;
w) span="$OPTARG";;
i) iteration="$OPTARG";;
x) skip="$OPTARG";;
d) redo="$OPTARG";;
n) label="n";;
p) path="$OPTARG";;
?) echo 'Unknown parameters.' && params && exit 1;;
esac
done

if [ -z "$gene" ] || [ -z "$species" ] || [ -z "$rtype" ] || [ -z "$bulks" ] || [ -z "$refline" ]
then
   echo -e '\nSome or all of the required parameters are empty.';
   params && exit 1
fi

if [ "$skip" = "help" ]
then
   echo -e '\nThis script overwrites files when rerun. This can be manually manipulated.'
   echo -e '\t-x option lets you skips some steps. This is useful when resuming analysis after an error.'
   echo -e '\t-x recognizes characters from e to m. To skip step 7 and 9, use "-x gi".'
   echo -e '\tUse these commands only when you know what you are doing.'
   echo -e '\t\te\t Setup snpEff custom annotation database'
   echo -e '\t\te\t Alignment of fastq file and BAM generation'
   echo -e '\t\tf\t BAM fixmate and sorting'
   echo -e '\t\tg\t BAM markduplicate, read group editing and indexing'
   echo -e '\t\th\t Variant calling'
   echo -e '\t\ti\t Calculate ratios for 2 bulks'
   echo -e '\t\tj\t Annotation of vcf files'
   echo -e '\t\tk\t Table generation' 
   echo -e '\t\tl\t Manhatten plot generation'
   echo -e '\t\tm\t CAPS marker identification' 
   echo -e '\t\tn\t Files cleanup'  
   exit 1
fi

if [ "$redo" = "help" ]
then
   echo -e '\nThis script recognizes and skips some repetative step based on presence or absence of certain files.'
   echo -e '\t-d option lets you force redo some steps. This is useful when rerunning an analysis or when corrupt outfile prevents the script from rerunning the step.'
   echo -e '\t-d recognizes characters from b to d. To redo step 3 and 4, use "-d cd".'
   echo -e '\t\tb\t Setup reference genome.'
   echo -e '\t\tc\t Setup known snps.'  
   echo -e '\t\td\t Setup custom annotation database.'
   exit 1
fi

########## Wrapper function for calling individual sub-analysis functions ##########

main() {

filesetup

export refpath="$refdir/$species/$refline"
[[ $redo =~ "b" || ! -f "$refpath/$refline.chrs.dict" ]] && refgenome

[[ $redo =~ "c" || ! -z "$snpdl" ]] && knownsnp

[[ $redo =~ "d" || `grep -q "$refline.genome : $species" $snpeffdir/snpEff.config` ]] && annotationbuild

[[ ! $skip =~ "e" ]] && module load parallel && parallel --jobs 2 alignment ::: "${bulks[@]}" ::: "${lanes[@]}"

[[ ! $skip =~ "f" ]] && module load parallel && parallel --jobs 2 bamprep ::: "${bulks[@]}" ::: "${lanes[@]}"

[[ ! $skip =~ "g" ]] && module load parallel && parallel --jobs 2 --tmpdir ./tmp bamprep2 ::: "${bulks[@]}"

[[ ! $skip =~ "h" ]] && variants

[[ ! $skip =~ "i" ]] && getratio "${bulks[0]}" "${bulks[1]}"

[[ ! $skip =~ "j" && $label = "y" ]] && annotation

[[ ! $skip =~ "k" && $label = "y" ]] && tables

[[ ! $skip =~ "l" ]] && plots

[[ ! $skip =~ "m" ]] && markers

[[ ! $skip =~ "n" ]] && cleanup

echo 'All done' && exit 1

}

########## Preparing input files ##########

filesetup() {

# Create required folders
[[ ! -d "fastq" && `ls | grep -E 'fastq.gz$|fastq$|fq.gz$|fq$' | wc -l` -eq 0 ]] && echo 'Error: No fastq files found.' && exit 1

[[ ! -d "fastq" ]] && mkdir fastq && mv *.fq *.fastq *.fq.gz *.fastq.gz fastq

[[ ! -d "$output" ]] && mkdir $output $output.archive

# Set temporary folder
export TMPDIR=./tmp

[[ ! -d "tmp" ]] && mkdir tmp
export _JAVA_OPTIONS=-Djava.io.tmpdir=tmp

echo 'Finished setting up the environment.'

}

########## Working with reference genome ##########

refgenome() {

module load samtools gatk bwa
[[ ! `type -t gatk` || ! `type -t samtools` || ! `type -t bwa` ]] && echo 'Error: Some modules not found.' && exit 1

# Create folder if necesary
mkdir -p $refpath

# Download and extract reference genome in ./refs/refline folder, skip if already done.
if [[ -z "$refdl" ]]
then
    echo -e '\n Please input reference genome link.'
	params && exit 1
else
    curl -o $refpath/$refline.fa.gz $refdl
    gzip -d $refpath/$refline.fa.gz
fi

# Separate chromosomes
awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' $refpath/$refline.fa > $refpath/$refline.chrs.fa

# Remove original fasta
rm $refpath/$refline.fa

# Create chromosome list from reference genome
samtools faidx $refpath/$refline.chrs.fa

# Indexing reference genome
bwa index $refpath/$refline.chrs.fa

# Create dictionary
gatk CreateSequenceDictionary -R $refpath/$refline.chrs.fa -O $refpath/$refline.chrs.dict

module purge

}

########## Known SNPs list preparation ##########

knownsnp () {

# Prepare known snp list
[[ -z "$snpdl" ]] && curl -o $refpath/knownSNPs.vcf.gz $snpdl && gzip -d $refpath/knownSNPs.vcf.gz

}

########## Alignment and BAM generation ##########

alignment() {

module load bwa samtools
[[ ! `type -t samtools` || ! `type -t bwa` || ! `type -t parallel` ]] && echo 'Error: Some modules not found.' && exit 1


nfiles=`ls fastq | grep ".*$1.*$2.*R[1-2].*" | wc -l`

if [ $nfiles -lt 1 -o $nfiles -gt 2 ]
then
	echo 'Skipping '$1.$2' because there are '$nfiles' fastq files. Expected 1 or 2 files.'
	[[ $nfiles -gt 2 ]] && echo 'Concatenate files before running the script.'
	return
fi
	
pool=( fastq/*$1.*$2*R[1-2].* )

# Alignment
bwa mem -t 8 -M $refpath/$refline.chrs.fa ${pool[*]} > $output/$gene.$1.$2.sam

# Binary file generation
samtools view -bSh $output/$gene.$1.$2.sam > $output/$gene.$1.$2.bam

module purge

}
export -f alignment

########## BAM manipulation ##########

bamprep() {

module load samtools
[[ ! `type -t samtools` || ! `type -t parallel` ]] && echo 'Error: Some modules not found.' && exit 1

[[ ! -f $output/$gene.$1.$2.bam ]] && echo 'Warning: Skipping '$1.$2' because bam file not found.' && return

#Big size SAMs can be deleted now
rm $output/$gene.$1.$2.sam

# For paired end only
[[ `samtools view -H $output/$gene.$1.$2.bam | grep "bwa.*$1.*$2.*R2.f"` ]] && samtools fixmate $output/$gene.$1.$2.bam $output/$gene.$1.$2.fix.bam

[[ -f $output/$gene.$1.$2.fix.bam ]] && mv $output/$gene.$1.$2.bam $output/$gene.$1.$2.nofix.bam && mv $output/$gene.$1.$2.fix.bam $output/$gene.$1.$2.bam

# Sort by coordinates
samtools sort -o $output/$gene.$1.$2.sort.bam $output/$gene.$1.$2.bam

module purge

}
export -f bamprep

bamprep2() {

module load samtools gatk
[[ ! `type -t gatk` && ! `type -t samtools` && ! `type -t parallel` ]] && echo 'Error: Some modules not found.' && exit 1

for i in ${lanes[@]}
do
	[[ ! -f $output/$gene.$1.$i.sort.bam ]] && echo 'Warning: Some lanes might be missing for '$1' bulk'
done

# Merge lanes
#files=( $output/$gene.$1.*.sort.bam )
#samtools merge $output/$gene.$1.sort.bam ${files[*]} #$(printf $output/$gene'.'$1'.%q.sort.bam ' "${lanes[@]}")

# Mark duplicates
gatk MarkDuplicates -I $output/$gene.$1.sort.bam -O $output/$gene.$1.md.bam -M $output/$gene.$1.matrics.txt -ASO coordinate
#SO coordina should do sorting too!! Check for this


# Add header
gatk AddOrReplaceReadGroups -I $output/$gene.$1.md.bam -O $output/$gene.$1.rg.bam -LB $gene.$1 -PL illumina -SM $gene.$1 -PU run1 -SO coordinate
## Adding RG in bwa can eliminate this step?

# Build index
gatk BuildBamIndex -I $output/$gene.$1.rg.bam

# Remove obsolete bams
#[[ -f $output/$gene.$1.rg.bai ]] && find $output -name *$gene.*$1.*bam | grep -v 'rg' | xargs rm

module purge

}
export -f bamprep2

########## Variant calling ##########

variants() {

module load gatk parallel
[[ ! `type -t gatk` || ! `type -t parallel` ]] && echo 'Error: Some modules not found.' && exit 1

for i in ${bulks[@]}
do
	[[ ! -f $output/$gene.$i.rg.bam ]] && echo 'Error: BAM file for '$i' bulk not found.' && exit 1
done

# GEtting list of chromosomes for scatter-gather
chrlist=( `grep 'SN:CM' $refpath/$refline.chrs.dict | cut -f2 | cut -d ":" -f2` )

# Varaint calling by chromosome
parallel --jobs 4 gatk HaplotypeCaller -R $refpath/$refline.chrs.fa $(printf ' -I '$output/$gene'.%q.rg.bam' "${bulks[@]}") -O $output/$gene.{}.hc.vcf -L {} ::: "${chrlist[@]}"

# Consolidate chroomsomes
gatk GatherVcfs $(printf ' -I '$output/$gene'.%q.hc.vcf' "${chrlist[@]}") -O $output/$gene.hc.vcf

module purge

}

########## Select variants ##########

getratio() {

module load gatk
[[ ! `type -t gatk` ]] && echo 'Error: Some modules not found.' && exit 1

#Select 2 bulks for QTL analysis
[[ ${#bulks[@]} > 2 ]] && mv $output/$gene.hc.vcf $output/$gene.hc.all.vcf && gatk SelectVariants -R $refpath/$refline.chrs.fa -V $output/$gene.hc.all.vcf -O $output/$gene.hc.vcf -sn $gene.${bulks[0]} -sn $gene.${bulks[1]}

#Change to tabular faromat
gatk VariantsToTable -R $refpath/$refline.chrs.fa -V $output/$gene.hc.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -O $output/$gene.table

# Calculating ΔSNP-index
printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" "chr" "pos"  "ref" "alt" "${bulks[0]}.ref" "${bulks[0]}.alt" "${bulks[0]}.ratio" "${bulks[1]}.ref" "${bulks[1]}.alt" "${bulks[1]}.ratio" "delratio\n" > $output/$gene.ratio.txt
cat $output/$gene.table | awk 'BEGIN{OFS="\t"} $2~/^[0-9X]*$/ && $8 !~ /NA/ && $12 !~ /NA/ && $3 !~ /\./ && $4 !~ /\./ && $7 >= 5 && $11 >= 5 {split($6,a,","); split($10,b,","); if ( a[1]+b[1] > 0 ) print $1, $2, $3, $4, a[1], a[2], a[2]/$7, b[1], b[2], b[2]/$11, a[2]/$7 - b[2]/$11 }' >> $output/$gene.ratio.txt

module purge

}

########## Build custom annotation database ##########

annotationbuild() {

[[ -z "$refdl" || -z "$annotdl" ]] && echo 'Error: Download links for reference genome not found. Skipping annotation building.' && params && label="n" && return

module load snpeff
[[ ! `type -t snpEff` ]] && echo 'snpEff not found. Forcing nolabel mode.' && label="n" && return

# Download annotation file
mkdir -p $snpeffdir/data/$refline
curl -o $snpeffdir/data/$refline/genes.gff.gz $annot

# Download reference genome
curl -o $snpeffdir/data/genomes/$refline.fa.gz $refdl

# Add genome to database
echo -e "\n# $species database\n$refline.genome : $species" >> $snpeffdir/snpEff.config

# Build custom database
snpEff build -c $snpeffdir/snpEff.config -gff3 -v $refline

module purge

}

########## Annotation ##########

annotation() {

[[ ! -f $output/$gene.hc.vcf ]] && echo 'Error: Input VCF file not found for annotation' && exit 1

[[ ! -f $snpeffdir/snpEff.config ]] && echo 'Warning: snpEff config file not found. Forcing nolabel mode.' && label="n" &&return

module load snpeff
[[ ! `type -t snpEff` ]] && echo 'Warning: snpEff not found. Forcing nolabel mode.' && label="n" && return

snpEff -c $snpeffdir/snpEff.config $refline -s $output/snpEff_summary.html $output/$gene.hc.vcf > $output/$gene.se.vcf

module purge

}

########## Finding candidates ##########

tables() {

[[ ! -f $output/$gene.se.vcf ]] && echo 'Warning: Annotated VCF not found. Forcing nolabel mode' && label=n && return

# Initializing headers for tables
printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" "chr" "pos" "ref" "alt" "mutation" "gene" "SNP.loc" "SNP.nt" "SNP.protein" "${bulks[0]}.ref" "${bulks[0]}.alt" "${bulks[1]}.ref" "${bulks[1]}.alt" | tee $output/$gene.all_SNPs.txt $output/$gene.seg_SNPs.txt $output/$gene.orf_SNPs.txt

# Go back to snpEff output format
awk 'FNR==NR{a[$1$2];next};($1$2 in a)' $output/$gene.ratio.txt $output/$gene.se.vcf  | cut -f3,7 --complement | awk '$3!~/\./ && $4!~/\./' > $output/$gene.snp.tmp

# Ordering SNPs by location
sort -k1,1 -k2,2n $output/$gene.snp.tmp > $output/$gene.snpsort.tmp

# Formatting the list for R
awk 'BEGIN{OFS="\t"} {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("n.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' $output/$gene.snpsort.tmp >> $output/$gene.all_SNPs.txt

# Finding all candidates for marking markers
case $rtype in
r) cat $output/$gene.snpsort.tmp | awk 'BEGIN{OFS="\t"} $8~/^1\/1/ && ($9~/^0\/1/ || $9~/^1\/0/ || $9~/^0\/0/) && $2~/^[0-9X]*$/ {print $0}' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2]) print $0}' > $output/$gene.segsnp.tmp;;
d) cat $output/$gene.snpsort.tmp | awk 'BEGIN{OFS="\t"} ($8~/^0\/1/ || $8~/^1\/0/ || $8~/^1\/1/) && $9~/^0\/0/ && $2~/^[0-9X]*$/ {print $0}' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2]) print $0}' > $output/$gene.segsnp.tmp;;
c) cat $output/$gene.snpsort.tmp | awk 'BEGIN{OFS="\t"} $8~/^1\/1/ && $9~/^0\/0/ && $2~/^[0-9X]*$/ {print $0}' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2]) print $0}' > $output/$gene.segsnp.tmp;;
?) echo 'Only dominant and recessive genes can be processed.' && params && return;;
esac

# Remove known SNPs
[[ -f $refpath/knownSNPs.vcf ]] && awk 'FNR==NR{a[$1$2];next};!($1$2 in a)' $refpath/knownSNPs.vcf $output/$gene.segsnp.tmp > $output/$gene.segsnpnew.tmp && mv $output/$gene.segsnpnew.tmp segsnp.tmp

# Formatting all candidates list
awk 'BEGIN{OFS="\t"} {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("n.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' $output/$gene.segsnp.tmp | awk '(($10+$11)>4) && (($12+$13)>4)' >> $output/$gene.seg_SNPs.txt

# Finding genic candidates
effects=( "bidirectional_gene_fusion" "coding_sequence_variant" "disruptive_inframe_deletion" "disruptive_inframe_insertion" "duplication" "duplication" "exon_loss_variant" "exon_loss_variant" "exon_variant" "feature_ablation" "frameshift_variant" "gene_fusion" "inframe_deletion" "inframe_insertion" "initiator_codon_variant" "inversion" "miRNA" "missense_variant" "rearranged_at_DNA_level" "splice_acceptor_variant" "splice_donor_variant" "splice_region_variant" "start_lost" "start_retained" "stop_gained" "stop_lost" )

grep ${effects[@]/#/-e } $output/$gene.seg_SNPs.txt >> $output/$gene.orf_SNPs.txt

}

########## Running R code ##########

plots() {

module load R
[[ ! `type -t Rscript` ]] && echo 'Error: R not found. No plots generated.' && return

[[ ! -f $output/$gene.ratio.txt ]] && echo 'Error: Ratio list not found. No plots generated.' && return
[[ $label = "y" && ! -f $output/$gene.orf_SNPs.txt ]]  && echo 'Warning: Candidate list not found. Forcing nolabel mode.' && label=n

awk -v thr=$thres 'BEGIN{OFS="\t"} {if (NR==1 || $11 >= thr || -$11 >= thr ) print $1, $2, $7, $10, $11}' $output/$gene.ratio.txt > $output/$gene.ratio.tmp

Rscript $path/NGS-BSA-plot.R "$output" "$gene" "$thres" "$span" "$iteration" "$label"

module purge

}

########## Finding markers ##########

markers() {

[[ ! -f $output/$gene.seg_SNPs.txt ]] && echo 'Error: Segregant SNP list not found. Markers are not developed.' && return
[[ ! -f $output/$gene.seg_SNPs.txt ]] && echo 'Error: Segregant SNP list not found. Markers are not developed.' && return
module load samtools
[[ ! `type -t samtools` ]] && echo 'Error: Samtools not found. Markers are not developed.' && return

SIZE=10; SIZE2=75

RE=`awk '{ORS="|"}{print $2}' REdatabase`

printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" "chr" "pos"  "ref" "alt" "ref.seq" "alt.seq" "ref.RE" "alt.RE" "misc" > $output/$gene.markers.txt

{
	read
	while IFS=$'\t' read -r CHR POS REF ALT REST
	do
		[[ ${#REF} > 36 || ${#ALT} > 36 ]] && COMMENT="Possible_PCR_marker" || COMMENT=""
	
		OLD=`samtools faidx $refpath/$refline.chrs.fa $CHR:$(($POS-$SIZE2))-$(($POS+${#REF}-1+$SIZE2)) | grep -v '>' | awk '{ print toupper($0) }'`
		NEW=${OLD:0:$SIZE2}$ALT${OLD:$SIZE2+${#REF}:$SIZE2}
		
		if [[ `grep -Eo $RE <<< $NEW` || `grep -Eo $RE <<< $OLD` ]]
		then
			REO=()
			REN=()		
		
			while IFS=$'\t' read -r NAME PAT
			do
				PATF=`awk '{ $0=toupper($0); gsub("W","[AT]"); gsub("S","[CG]"); gsub("M","[AC]"); gsub("K","[GT]"); gsub("R","[AG]"); gsub("Y","[CT]"); gsub("B","[CGT]"); gsub("D","[AGT]"); gsub("H","[ACT]"); gsub("V","[ACG]"); gsub("N","."); gsub("Z",""); print }' <<< "$PAT"`
				PATR=`echo "$PATF" | tr "[ATGCUWSMKRYBDHV]" "]TACGAWSKMYRVHDB[" | rev`;
				[[ `grep -Eo "$PAT|$PATR" <<< $OLD | wc -l` -ge 1 && `grep -Eo "$PAT|$PATR" <<< $NEW | wc -l` -eq 0 ]] && REO+=($NAME)
				[[ `grep -Eo "$PAT|$PATR" <<< $OLD | wc -l` -eq 0 && `grep -Eo "$PAT|$PATR" <<< $NEW | wc -l` -ge 1 ]] && REN+=($NAME)
			done < $path/REdatabase
	
			OLD=`samtools faidx ../reference/Capsicum_annuum/CM334/CM334.chrs.fa $CHR:$(($POS-$SIZE))-$(($POS+${#REF}-1+$SIZE)) | grep -v '>'`
			NEW=${OLD:0:$SIZE}$ALT${OLD:$SIZE+${#REF}:$SIZE}
	
			printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" $CHR $POS $REF $ALT $REST $OLD $NEW $(IFS=','; echo "${REO[*]}"$"\t""${REN[*]/%/*}") $COMMENT >> $output/$gene.markers.txt
		fi
	done 
} < $output/$gene.seg_SNPs.txt

module purge

}

########## Cleanup ##########

cleanup() {

rm -r tmp
rm -r $output/*.tmp
mv $output/*.table $output.archive
mv $output/*.bam $output.bams
mv $output/*.vcf $output.vcf

}

main
