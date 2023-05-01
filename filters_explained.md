## Class: FilterVCF

This is the code I used to filter SNPs down after variant calling (before final mapping).   
Written by N. Schlenk unless otherwise stated.  
Scripts can be found in **toolbox** under class **FilterVCF**.

The **FilterVCF** class has 3 arguments:
```{}
out_dir          # the directory automatically assigned for every output file (does not apply to inputs, default wd)
nF2s             # number of F2s in the sample (default 238)
scaffold_index   # index for defining characters in the scaffold titles (this will only work for single digits or letters, default 12)
```

To run any of the following scripts, use **FilterVCF(out_dir = "", nF2s = "", scaffold_index = "")** followed by **.ModuleName(module args)**.

### Module: MQADPAR
Filter each SNP for MQ sccore, distribution, and parent genotypes. Written by C. Wessinger and J. Kelly, modified by N. Schlenk.

```{}
infile_path   # this is the VCF to be filtered 
outfile_path  # this is the name of the output VCF (do not include out_dir if one is included above)
minIndiv      # minimum number of individuals carrying SNP (default 50)
minMQ         # minimum MQ score (default 30)
```

### Module: BestSNP
This filter finds the best SNP per 150bp and filters for minor allele distribution. Written by J. Kelly, modified by N. Schlenk.

```{}
infile_path   # this is the VCF to be filtered
outfile_path  # this is the name of the output VCF (do not include out_dir if one is included above)
minMinor      # minimum number of individuals carrying minor allele
```

### Module: ParentDepth
Filter each SNP based on allele depth of parents (assumed to be the last 2 columns). *I did not use this filter in my project, but I'll keep it here in case it is useful.*
```{}
infile_path   # this is the VCF to be filtered
outfile_path  # this is the name of the output VCF (do not include out_dir if one is included above)
p1a1min       # minimum number of p1 alleles in the first parent (p1), default 100
p2a1min       # minimum number of p1 alleles in the second parent (p2), default 0
p1a2min       # minimum number of p2 alleles in the first parent (p1), default 100
p2a2min       # minimum number of p2 alleles in the second parent (p2), default 0
```

### Module: FilterHWE
Filter each SNP based on genotype frequencies (or counts). This module doesn't actually use HWE but can be helpful for narrowing down HWE-friendly SNPs.
```{}
infile_path   # this is the VCF to be filtered
outfile_path  # this is the name of the output VCF (do not include out_dir if one is included above)
aa_min        # minimum frequency (or count) of '0/0' genotypes
ab_min        # minimum frequency (or count) of '0/1' genotypes
bb_min        # minimum frequency (or count) of '1/1' genotypes
count         # boolean (true/false), determines if filter will apply to genotype counts per SNP (true) or genotype frequencies per SNP (false, default)
```
