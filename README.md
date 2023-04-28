## Project Overview

**Title:** The _Penstemon kunthii_ draft genome: Integrating a genetic map with assembled sequence data  
**Purpose:** To verify the existing _Penstemon kunthii_ genome assembly using a linkage map derived from a cross between _P. kunthii_ and _P. amphorellae_. This was my master's thesis project at the University of Kansas.  
**Summary:** We generated MSG sequence reads from an F2 population, which were then used to make a genetic map for each of the 8 _Penstemon_ chromosomes. Then, I investigated other genomic elements including GC content and repeat element distribution. This repository holds the code used to complete this project, as well as my notes/instructions.
  
**Required software:** anaconda (msg_ipyrad, NumPy), bwa, samtools, BCFtools, java, python3, Lep-MAP3  
**Required files:** picard.jar (attached), sequences/barcodes (multiplexed), pedigree

A description of my workflow can be found in **project_steps_explained.md** above.  
Custom scripts I used to filter, order, plot, and otherwise analyze the variant-called dataset are all included in **toolbox** (above).
