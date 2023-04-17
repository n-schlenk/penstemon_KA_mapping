##### Required software:
- bioconda (to install msg_ipyrad environment)
- anaconda
- bwa
- samtools
- bcftools
- java
- python3
- NumPy
- SciPy
- Lep-MAP3 (download from SourceForge)

##### Required files:
- picard.jar (attached)

## Demultiplex and Filter
|||
|-----|-----|
|Starts with:|FASTQ (seq output), TXT (barcodes)|
|Ends with:|FASTQ (concatenated by ID), TXT (stats)|
##### Workflow:
First, you must load the msg_ipyrad (bioconda) environment. Use this to create parameters file.
```
ipyrad -n <title>
```
Then, submit the job (recommended sbatch to cluster)
```
ipyrad -p <parameter file> -s 12 -c 1 -f
```
##### Flags and parameters:
```
ipyrad -n <title>                       # creates new parameter file, which must be manually edited
ipyrad -s 12                            # only complete steps 1 and 2 of ipyrad workflow (demultiplex, filter)
ipyrad -c 1                             # print results (if you sbatch, these will end up in a log file)
ipyrad -f                               # forces ipyrad to run steps, even if they have already been completed (applies when some data is missing)                         
```

## Align to Reference Genome
|||
|-----|-----|
|Starts with:|FASTQ (one per ID), FASTA (ref genome)|
|Ends with:|BAM (one per ID)|
##### Workflow:
```
bwa index | bwa mem | java AddOrReplaceReadGroups | samtools view -b | samtools sort > <final BAM file>
```
##### Flags and parameters:
```
samtools view -b                        # translates the SAM to BAM (binary)
```
## Call Variants
|||
|-----|-----|
|Starts with:|BAM (one per ID)|
|Ends with:|VCF (single)|
##### Workflow:
```
bcftools mpileup | bcftools call > <VCF>
```
##### Flags and parameters:
```
bcftools mpileup -Ou                     # output to uncompressed BCF
bcftools mpileup -I                      # exclude indels (only include SNPs)
bcftools mpileup -a FORMAT/AD            # outputs allelic depth in addition to the defaults (genotype and phred-scaled genotypic likelihood) 
bcftools mpileup --max-depth=100         # include a maximum of 100 reads per sample per SNP
```
```
bcftools call -v                         # output variant sites only 
bcftools call -m                         # use alterate model for multiallelic and rare-variant calling (recommended)
bcftools -O v                            # output to uncompressed VCF
```

## Additional Filtering
|||
|-----|-----|
|Starts with:|VCF, pedigree|
|Ends with:|better VCF|

java -cp <filepath to Lep-MAP3 bin> is required before ParentCall2 and Filtering2 because they require Lep-MAP3 information
##### Filters (complete in order):
```
python filter_mq_aldepth_parents.py                     # removes markers with low MQ score, low call counts, more than one alt allele, and those which display noninformative/unexpected parental genotypes (parental genotypes must be the last 2 columns)
python filter_refparent.py                              # removes markers where the parental genotype of the genome species is not homozygous for the ref allele
python filter_bestsnp.py | filter_extractsnp.py         # filters VCF for best SNPs per radtag (assumes that you have already filtered for quality)
head -n37 unfiltered.vcf > pre_pc.vcf                   # add original headers back to VCF before using ParentCall2
cat filtered.vcf >> pre_pc.vcf
java -cp <LM filepath> ParentCall2 data=pedigree.txt vcfFile=filtered.vcf > pc.vcf

```
##### Parameters:
`````
filter_mq_aldepth_parents.py                            # minIndiv = 50 (minimum number of called genotypes at the SNP)      
filter_mq_aldepth_parents.py                            # minMQ = 30 (minimum MQ score for the SNP)
filter_bestsnp.py                                       # MinMinor = 1 (number of individuals that will have the minor allele)
ParentCall2 removeNonInfomative=1                       # removes markers that are monomorphic or homozygous for both parents
`````

## Map
|||
|-----|-----|
|Starts with:|filtered and parentcalled VCF|
|Ends with:|genetic map coordinates for each scaffold|
##### Workflow:
```
sh separate_scaffold filtered_parentcalled.vcf          generates VCF and map file for each scaffold
sh ordermarkers.sh                                      applies Lep-MAP3 OrderMarkers2 on each scaffold
python parentcall_reformat.py                           create new reference VCF from parentcalled VCF file
python get_coords.py                                    generates coordinates for output of OrderMarkers2 (per scaffold)
```
##### Additional Arguments:
```
OrderMarkers2 sexAveraged=1                             cM values for markers are averaged between sexes
OrderMarkers2 improveOrder=0                            calculate cM values with the orders in which SNPs are aligned
```

### Pre-Graphing
|||
|-----|-----|
|Starts with:|map coordinates|
|Ends with:|data for plotting|
##### Workflow:
```
python find_gc_per_chr.py                                 input fastas, output txt files with GC content per <bin> 
python find_repeats_per_chr.py                            input stockholm file (RepeatModeler) and target repeat element family, output repeat element count per <bin>

```

## Demultiplex 
*This process is specific to the way that we have been using ipyrad at KU, but the steps could be applied to similar projects and environments. Parameters were provided by Dr. Carrie Wessinger and Dr. Lena Hileman.*   

To start, we need a **raw fastq** file and a **barcode** file. 

**Raw fastq** files should have the **fastq.gz** extension. Illumina explains these files really well on their website: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html. These files are often named in an unusual way and might look a little confusing. An example is shown below.

```{}
LCH-AseI-P1_S1_R1_001.fastq.gz
```
  
**Barcode** files should be in **.txt** format and contain only 2 columns separated by a single tab. The first column should contain the unique specimen ID and the second should contain its barcode sequence. A 5-specimen example is shown below.

```{}
AB001   CGTGG
AB002   GTTGG
AB003   CCATG
AB004   GGTAA
AB005   CCTCA
```
  
Once these files are ready to go (in a community cluster, they should reside in a **temp** or **scratch** folder), we can start to demultiplex them. I recommend doing this in an interactive job.

```{}
srun --time=4:00:00 --ntasks=1 --nodes=1 --partition=sixhour --pty /bin/bash -l
```
  
Then, we are going to load **conda** so that we may proceed in a **conda** environment called **MSG ipyrad**. This will take some time to complete.
  
```{}
module load conda
conda activate msg_ipyrad
ipyrad -n demultiplex-plate1
```
   Note that in this step, we also used ipyrad to generate a **parameters** file. Above, it is named "demultiplex-plate1", but it can be anything. This new **parameters** file needs to be edited. I used **nano** in the command line.

```{}
nano parameters-demultiplex-plate1
```
  
There are several parameters that need to be changed. Here is what our final **parameters** file looks like. 
```{}
------- ipyrad params file (v.0.9.81)-------------------------------------------
demultiplex-plate1                  ## [0] [assembly_name]: Assembly name. Used to name output directories...
/temp/30day/hileman/n752s925/KA.MSG ## [1] [project_dir]: Project dir (made in curdir if not present)
LCH-AseI-P1_S1_R1_001.fastq.gz      ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
apr2021.AseI.plate1.barcodes.txt    ## [3] [barcodes_path]: Location of barcodes file
                                    ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                              ## [5] [assembly_method]: Assembly method (denovo, reference)
                                    ## [6] [reference_sequence]: Location of reference sequence file
gbs                                 ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TAAT,                               ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                                   ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                                  ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default...
6                                   ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                                   ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                               ## [13] [maxdepth]: Max cluster depth within samples
0.85                                ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                                   ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                                   ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
75                                  ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                                   ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
5                                   ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
8                                   ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
4                                   ## [21] [min_samples_locus]: Min # samples per locus for output
20                                  ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                                   ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                                 ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                          ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                          ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                             ## [27] [output_formats]: Output formats (see docs)
                                    ## [28] [pop_assign_file]: Path to population assignment file
                                    ## [29] [reference_as_filter]: Reads mapped to this reference are removed...
```
It is also important to remember that, at this step, we copy the name of the **raw fastq** file and **barcode** file exactly.  


  
 Next, we can run the job. I did this by submitting it to our department's partition in our community cluster using a shell script. I generated the shell script using **nano**, copied and pasted the following code, and then used **sbatch** to submit it. If this is not being completed on a community cluster, only the last line should be used.
   
```{}
#!/bin/bash
#SBATCH --job-name=ipyrad_demultiplex    # Job name
#SBATCH --partition=eeb                  # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=n752s925@ku.edu      # Where to send mail
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=8                # Number of CPU cores per task
#SBATCH --mem=20gb                       # Job memory request
#SBATCH --time=0-24:00:00                # Time limit days-hrs:min:sec
#SBATCH --output=out_%j.log              # Standard output and error log

cd $SLURM_SUBMIT_DIR

module load anaconda
conda activate msg_ipyrad

ipyrad -p params-demultiplex-plate1.txt -s 12 -c 1 -f
```
This code utilizes the following flags:
```{}
ipyrad ... -s 12                       # only complete steps 1 and 2 of ipyrad workflow (demultiplex, filter)
ipyrad ... -c 1                        # print results (these will end up in a log file when using sbatch)
ipyrad ... -f                          # forces ipyrad to run steps, even if they have already been completed (applies when some data is missing)      

```
When the code is done, there are 2 new directories:
  
```{}
demultiplex-plate1_edits
demultiplex-plate1_fastqs
```

In the **edits** folder, there are all **trimmed R1 fastq.gz** files and a **stat** file, both of which are important. We keep everything from the **edits** folder.  In the **fastqs** folder, we only keep the **stat** file. the **untrimmed/original quality fastq** files can be deleted.  

It is recommended to move **everything** from the **edits** folder and the **stat** file from the **fastqs** folder to long-term storage (we use ResFS).  




## Align to Reference Genome

Once all of the **demultiplexed** fastq files are in the same directory (here called "KA_fastqs"), we may proceed with aligning them. This step requires the following modules: bwa, java, picard, samtools. We also need the **reference genome** in fasta format (here called "PGA_assembly_trimmed.fasta").
The first step is to index the **reference genome**.

```{}
bwa index -p kunthii PGA_assembly_trimmed.fasta
```

Then, we loop through each of the **fastq** files to produce an alignment file in the sequence alignment map (**SAM**) format. 


In this example, I allowed the **sample ID** to be defined by the first 5 characters of the **fastq** filepath (see line 3). This same **sample ID** is also used to define a **read group**, which we will add to an intermediate **SAM** file before we sort and translate it to generate our output in the binary alignment map (**BAM**) format. 


```{}
for i in $(ls KA_fastqs/)
do
    ID=$(echo $i | grep -Po '^.{5}')
    bwa mem kunthii KA_fastqs/$i > $ID'_bwa_output.sam'                 # the alignment step
    java -jar picard.jar AddOrReplaceReadGroups \                       # this adds read groups to the SAM file
        I=$ID'_bwa_output.sam'\
        O=$ID'_withRG.sam' \
        RGLB=$ID \
        RGPL=ILLUMINA \
        RGPU=barcode \
        RGSM=$ID
    samtools view -b $ID'_withRG.sam' | samtools sort > $ID'.bam'       # translate to BAM (-b) and sort by leftmost coordinates
done
```

The addition of a **read group** is not always necessary. In this case, we arbitrarily named each **read group** after the **sample ID**, which corrected a downstream error generated by a lack of **read groups**. It is possible that this step (the code following "java" and preceeding "samtools") can be skipped.


## Call Variants

Now that we have alignment information (**BAM** files) for each **fastq**, we can call variants (otherwise known as "calling SNPs"). This step requires **bcftools**. We also need a list (**.txt**) of all target **BAM** files. Mine are kept in a directory called "bam_files", so I generated my list (bam_list.txt) like this:

```{}
ls bam_files > bam_list.txt
```

Now, we can perform variant calling. There are 2 steps, which will be described seperately but I recommend piping them together. The first step is **mpileup**, which outputs an intermediate binary variant call format (**BCF**) file. In this step, **fastq** files are combined based on **read group** (which in this case just organizes them by **sample ID**).

```{}
bcftools mpileup -Ou -I -a FORMAT/AD --max-depth=100 -f PGA_assembly_trimmed.fasta -b bam_list.txt
```

The next step applies an algorithm to call variants on the **BCF** file. With these flags, the output is a variant call format (**VCF**) file (which I have named "unfiltered.vcf").

```{}
bcftools call -v -m -O v > unfiltered.vcf
```

The flags used in these 2 steps are:
```{}
bcftools mpileup ... -Ou                     # output to uncompressed BCF
bcftools mpileup ... -I                      # exclude indels (only include SNPs)
bcftools mpileup ... -a FORMAT/AD            # outputs allelic depth in addition to the defaults (genotype and phred-scaled genotypic likelihood) 
bcftools mpileup ... --max-depth=100         # include a maximum of 100 reads per sample per SNP
bcftools call ... -v                         # output variant sites only 
bcftools call ... -m                         # use alterate model for multiallelic and rare-variant calling (recommended)
bcftools ... -O v                            # output to uncompressed VCF
```

## Generate Linkage Map

*Note: Before making a linkage map, it might be necessary to apply filters.*
We used **Lep-MAP3** for linkage mapping. There are lots of modules in this software, but we won't use all of them.
To order markers, we use the OrderMarkers2 function. This requires a VCF generated from ParentCall2, so our first step is to "parentcall" the VCF we plan to use for linkage mapping. This eliminates markers based on genotype likelihood. The algorithm for this function can be found in the software documentation.

```{}
ParentCall2 data=pedigree.txt vcfFile=filtered.vcf > parentcalled.vcf
```

After we have filtered down and "parentcalled" the VCF, we can order the makers. We can do this with or without a genome. 
If we want to calculate genetic distance using the order of SNPs in the genome, we need to generate a "fake" map for OrderMarkers2. I found it easier to do this per-chromosome, so I wrote a script (Order.OrderMarkersGenome) which separates the parentcalled VCF into scaffolds and makes the maps for OrderMarkers2, then calls OrderMarkers2. An example line for OrderMarkers2 is shown below. This is executed as part of the OrderMarkersGenome function.

```{}
OrderMarkers2 evaluateOrder=map.txt data=parentcalled.vcf sexAveraged=1 improveOrder=0 > ordered.txt
```

The parameters used in this step are:
```{}
OrderMarkers2 ... evaluateOrder=<file>             # input map file (integers from 1 to the number of markers, corresponding to marker order)
OrderMarkers2 ... data=<file>                      # input parentcalled VCF
OrderMarkers2 ... sexAveraged=1                    # average genetic distances calculated for each parent
OrderMarkers2 ... improveOrder=0                   # do not rearrange markers
```

To generate genetic distances without using the ordering in the genome, it is typical to have Lep-MAP3 split the entire parentcalled VCF into linkage groups based on a LOD parameter. For the purposes of this project, I only had Lep-MAP3 order the SNPs within genome-determined linkage groups. I wrote a script (Order.OrderMarkersLepMAP) which splits the parentcalled VCF and runs OrderMarkers2 as below for each linkage group. The variable **i** defines the current linkage group and **n** determines the iteration. It is recommended to repeat OrderMarkers2 a few times (I did 3).

```{}
for i in chromosomes:
    OrderMarkers2 map=map.txt data=parentcalled.vcf sexAveraged=1 chromosome=i > ordered_i.txt
    most_recent_iteration = 'ordered_i.txt'
    for n in iterations:
        OrderMarkers2 evaluateOrder=most_recent_iteration data=parentcalled.vcf sexAveraged=1 chromosome=i > ordered_i.n.txt
        most_recent_iteration = 'ordered_i.n.txt'
```

The parameters used here are:
```{}
OrderMarkers2 ... map=<file>                       # input map file (assigns markers to linkage group)
OrderMarkers2 ... data=<file>                      # input parentcalled VCF
OrderMarkers2 ... sexAveraged=1                    # average genetic distances calculated for each parent
OrderMarkers2 ... chromosome=<n>                   # determine which linkage group to order
OrderMarkers2 ... evaluateOrder=<file>             # improve the order of the previous output of OrderMarkers2
```

