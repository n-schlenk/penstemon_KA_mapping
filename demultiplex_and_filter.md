The purpose of this document is to simplify the process of demultiplexing and assembling MSG reads.   
*This process is specific to the way that we have been using ipyrad at KU, but the steps could be applied to similar projects and environments. The code was generated by Dr. Carrie Wessinger and Dr. Lena Hileman.*   

For more detail on the **ipyrad** workflow and parameters, please see **ipyrad_parameters_explained**.

To start, you need a **raw fastq** file and a **barcode** file. Sometimes they aren't formatted in the correct way, so it's important to look at the data before proceeding.  

**Raw fastq** files should have the **fastq.gz** extension. Illumina explains these files really well on their website: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html. These files are often named in a funky way and might look a little confusing. An example is shown below.

```{}
LCH-AseI-P1_S1_R1_001.fastq.gz
```
  
**Barcode** files should be in **.txt** format and contain only 2 columns separated by tab. The first column should contain the unique specimen ID and the second should contain its barcode sequence. A 5-specimen example is shown below.

```{}
AB001   CGTGG
AB002   GTTGG
AB003   CCATG
AB004   GGTAA
AB005   CCTCA
```
  
Once you have these files in your **temp** or **scratch** folder, we can start to demultiplex them. Log onto the CRC and nagivate to your **temp** folder to ensure that your files are there. Then, start an interactive job.

```{}
srun --time=4:00:00 --ntasks=1 --nodes=1 --partition=sixhour --pty /bin/bash -l
```
  
Then, we are going to load **anaconda** so that we may proceed in a **conda** environment and utilize the **MSG ipyrad** tool. Running the code below will take some time to complete.
  
```{}
module load anaconda
conda activate msg_ipyrad
ipyrad -n demultiplex-plate1
```
   Note that in this step, we also used ipyrad to generate a **parameters** file. Above, it is named "demultiplex-plate1", but you can change this depending on what you're doing. This new **parameters** file needs to be edited. To edit, use the **nano** function. 

```{}
nano parameters-demultiplex-plate1
```
  
There are several parameters that need to be changed. See the default parameters below for comparison. 
```{}
------- ipyrad params file (v.0.9.81)-------------------------------------------
demultiplex-plate1             ## [0] [assembly_name]: Assembly name. Used to name output directories...
/temp/30day/hileman/n752s925/KA.MSG ## [1] [project_dir]: Project dir (made in curdir if not present)
LCH-AseI-P1_S1_R1_001.fastq.gz  ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
apr2021.AseI.plate1.barcodes.txt   ## [3] [barcodes_path]: Location of barcodes file
                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
gbs                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TAAT,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default...
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
75                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
4                              ## [21] [min_samples_locus]: Min # samples per locus for output
20                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                              ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                      ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
                               ## [29] [reference_as_filter]: Reads mapped to this reference are removed...
```
It is also important to remember that, at this step, you are copying the name of the **raw fastq** file and **barcode** file exactly.  
  
  Next, we will make a **shell script** file, which has the **.sh** extension. This is a script that we will later instruct the supercomputer to execute. If you ask it to demultiplex 2 plates at once, you will need to modify the **name** and **parameter** file of this **script**.  
  
```{}
nano ipyrad.sh
```
  
  Once you make the empty **shell script** file, paste the following information (all symbols included) into the document. Make sure to change the *where to send mail* option to include your KU email address. It is also important to copy the name of your **parameters** file exactly as it is.  

```{}
#!/bin/bash
#SBATCH --job-name=ipyrad_demultiplex    # Job name
#SBATCH --partition=eeb           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=n752s925@ku.edu        # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=8             # Number of CPU cores per task
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=0-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=out_%j.log   # Standard output and error log

cd $SLURM_SUBMIT_DIR

module load anaconda
conda activate msg_ipyrad

ipyrad -p params-demultiplex-plate1.txt -s 12 -c 1 -f
```

Now, submit the job to the cluster.
```{}
sbatch ipyrad.sh
```
  
  At this point, you can close your window and exit the cluster. You will get emails as your job progresses. When it is done, you will see 2 new directories:
  
```{}
demultiplex-plate1_edits
demultiplex-plate1_fastqs
```

In the **edits** folder, you will have all **trimmed R1 fastq.gz** files and a **stat** file, both of which are important. Keep everything from the **edits** folder.  In the **fastqs** folder, you will only need to keep the **stat** file. You can delete all of the **untrimmed/original quality fastq** files.  

Move **everything** from the **edits** folder and the **stat** file from the **fastqs** folder to ResFS.  
