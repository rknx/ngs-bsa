# ngs-bsa
Pipeline for the Next Gen Sequencing Bulked Sequence Analysis pipeline

### Instructions:
- Install the following programs: Java, Samtools, BWA, SnpEff, GATK and R.
- Install following packages in R: ggplot2, ggrepel, reshape2.
- Place fastq files (fastq.gz, fastq, fq.gz or fq) from both pools in root folder.
- Run script with appropriate parameters.
- If using SLURM, fill in the variables in batch file and run this script from batch file (see example batch file).

### Usage:
```
NGS_BSA.sh [-h <help>] -g gene -s species [-o output folder] -t resistance type -b {bulk codes} {[-l lane names]} [-r reference folder] -v reference line [-f reference genome link] [-a annotation link] [-e snpEff config folder] [-k known snps link] [-c cut-off threshold] [-w smoothing span] [-i iteration] [-x steps to skip] [-d steps to redo] [-n <nolabel>] [-p path to script folder]

where:  
        -g        Short name of gene under analysis (alphanumeric)  
        -s        Species name, (alphanumeric and "_")  
        -o        Output folder (alphanumeric and "_")  
        -t        Resistance type ("d: dominant", "r: recessive" or "c: codominant")  
        -b        List of names or codes for bulks or samples (alphanumeric and "_").  
                        Input resitant bulk followed by susceptible bulk.  
                        Only first two bulks will be used for plotting.  
        -l        Lane names if done on multiple lanes (alphanumeric and "_")  
        -r        Reference directory path (Created if not present)  
        -v        Reference variety (alphanumeric and "_")  
        -f        Reference genome fasta gzip download link  
        -a        Annotation file gzip link for reference genome  
        -e        snpEff config file path  
        -k        Known SNPs vcf gzip download link  
        -c        Cut-off threshold for manhattan plot (numeric, 0-1, default = 0.3)  
        -w        Spanning width for smoothing in Manhattan plot (numeric, 0-1, default = 0.2)  
        -i        Run iteration (numeric, default=NULL)  
        -x        Manually skip some steps. Run -x help for list of steps  
        -d        Force redo some steps. Run -d help for list of steps  
        -n        Remove gene label in plots (numeric, default=NULL)  
        -p        Path to script folder  
  ```

### Example:
        See the example batch file.

### Citation:
        Under review. To be listed soon.  
