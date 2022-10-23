#!/usr/bin/env python
# coding: utf-8

# # Task 2:  ChIP-seq data analysis
# ## Notebook
# 
# Comparative & Regulatory Genomics - [I0U29a] | Task 2 - ChIP-seq data analysis | Antoine Ruzette <b> r0829308 </b> | 19.12.2021
# 
# Based on the ChIP-seq jupyter notebook from Prof. Stein Aerts. 

# ### Setting up the environment

# In[3]:


mkdir -p /mnt/storage/$USER/jupyternotebooks/Assignment2
cd /mnt/storage/$USER/jupyternotebooks/Assignment2


# In[4]:


vdb-config -s /repository/user/cache-disabled=true


# ### Data loading

# ##### ChIP-seq data: Six1_myoblasts

# In[5]:


#data loading
fastq-dump --split-files SRR14711117

#delete it after the mapping done (should be about 30GB)


# ##### ChIP-seq data: Input_myoblasts

# In[12]:


fastq-dump --split-files SRR14711120


# In[22]:


ls -lt *.fastq


# ### Quality Check (QC)

# In[13]:


fastqc SRR14711117_1.fastq
fastqc SRR14711120_1.fastq
#insert screenshot for the three quality metrics


# In[21]:


ls -lt *.zip


# ### Aligment of the Six1 ChIP-seq reads from mouse primary myoblast cells

# To align the reads to the genome, we will use bowtie2. As reference genome, you can use the whole genome bowtie index at this location /mnt/storage/data/resources/mm10/mm10. Alternatively, you can use bwa-mem to align reads. Here we will align to mm10 mouse assembly. 

# In[11]:


#mapping the mouse genome
bowtie2 -x /mnt/storage/data/resources/mm10/mm10 SRR14711117_1.fastq -S ChIP_Six1.sam  


# In[20]:


ls -l *.sam


# Check a couple of lines in the SAM file

# In[25]:


head -500 ChIP_Six1.sam | tail -5


# Convert the SAM file to BAM format (binary)

# In[17]:


samtools view -S -b ChIP_Six1.sam > ChIP_Six1.bam


# How many reads are in the BAM file?

# In[30]:


samtools view -c ChIP_Six1.bam


# In[31]:


# List all SAM flags.
samtools flags


# In[32]:


# Get correct flag settings for samtools view (second column).
samtools flags UNMAP,SECONDARY


# In[33]:


# Number of mapped reads.
samtools view -c -F 260 ChIP_Six1.bam


# Now we are going to compress the .sam into .bam files and create the .bai index so we can use the BAM file in IGV

# In[34]:


samtools sort -O bam -o ChIP_Six1.sorted.bam ChIP_Six1.bam


# In[35]:


samtools index ChIP_Six1.sorted.bam


# In[43]:


ls ChIP_Six1*
#Note that we now have four files: sam > bam > sorted.bam > sorted.bam.bai


# ### Check your BAM file in IGV

# Copy the `ChIP_Six1.sorted.bam` and `ChIP_Six1.sorted.bam.bai` to your computer and open them with IGV. Go to the control genes (type the gene name in the search box and press "Go").

# ### Alignment of the control data ("input myoblastes")

# In[ ]:


bowtie2 -x /mnt/storage/data/resources/mm10/mm10 SRR14711120_1.fastq -S ChIP_input.sam


# In[ ]:


samtools view -S -b ChIP_input.sam > ChIP_input.bam


# How many reads are in the BAM file?

# In[1]:


samtools view -c ChIP_input.bam


# In[2]:


# Get correct flag settings for samtools view (second column).
samtools flags UNMAP,SECONDARY


# In[3]:


# Number of mapped reads.
samtools view -c -F 260 ChIP_input.bam


# Now we are going to create the .bai index so we can use the BAM file in IGV

# In[4]:


samtools sort -O bam -o ChIP_input.sorted.bam ChIP_input.bam


# In[5]:


samtools index ChIP_input.sorted.bam


# In[2]:


ls ChIP_input*
#go on IGV, insert screenshots - no need for the input file


# ### Genome-wide coverage plots

# We will generate a bigwig file that contains only the coverage of the reads, so not the individual reads. This bigWig file is a binary file format (simlar to BAM), and it can be accessed over the internet. This means that you can use the bigWig file locallly in IGV, but also as custom track in UCSC (provided that the bigWig is available from a web/ftp server, so you can specific the custom track by entering the URL).

# We're using the command bamCoverage from the package deeptools to create a bigwig file that can be used for visualization in IGV or UCSC Genome Browser. 
# 
# 
# Find the size of the whole mm10 mouse genome in UCSC Genome Browser. Replace the green value by the latter. 

# In[8]:


bamCoverage -b ChIP_Six1.sorted.bam --normalizeUsing RPGC --effectiveGenomeSize 2730855475 -o six1.bw


# Open that .bw file in IGV (=> green track)

# Finally, from the ChIP-atlas page of this data set http://chip-atlas.org/view?id=SRX285808, I have chosen "View on IGV - bigWig", which automatically opens the bigWig file in my local IGV (blue track)
# 

# ### Peak calling 
# 

# Next, we will determine genomic regions with an enrichment of reads, in the form of a "peak". 

# If you color the reads by strand (right click in IGV, color alignments by read strand), you will notice the shift between forward and reverse reads.

# You can also group reads by read strand - here is another example of a target gene, Slc16a10, located on the 10th chromosome.  Again we observe the shift in reads between positive and negative strand.
# 
# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# ### Peak calling with MACS2

# MACS implements a model to take advantage of the strand shift, to distinguish true, bona fide peaks, from artifacts.
# 
# Genome size for the reference genome mm10 was found in [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/). 

# In[4]:


macs2 --version


# In[ ]:


macs2 callpeak -t ChIP_Six1.sorted.bam -c ChIP_input.sorted.bam -n six1 -g 2730855475 -q 0.5
#replace genome size
#-c ChIP_input.sorted.bam
#maybe change the qvalue threshold


# In[7]:


ls *peaks*


# macs2 produces the following 4 files:
# 
# * `Six1_peaks.xls`: is a tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions. Information include position, length and height of detected peak etc.
# * `Six1_peaks.narrowPeak`: is BED6+4 format file which contains the peak locations together with peak summit, p-value and q-value. You can load it directly to IGV or UCSC genome browser.
# * `Six1_summits.bed`: is in BED format, which contains the peak summits locations for every peaks. The 5th column in this file is -log10p-value the same as NAME_peaks.bed. If you want to find the motifs at the binding sites, this file is recommended. The file can be loaded directly to UCSC genome browser. But remember to remove the beginning track line if you want to analyze it by other tools.
# * `Six1_model.r`: is an R script which you can use to produce a PDF image about the model based on your data. Load it to R by: `$ Rscript NAME_model.r` Then a pdf file NAME_model.pdf will be generated in your current directory. Note, R is required to draw this figure.

# How many peaks were called?

# In[9]:


cat six_1_peaks.narrowPeak | wc -l


# Copy the `Six1_peaks.narrowPeak` file to your computer, and open it in IGV.

# In[7]:


sort -k 9 six_1_peaks.narrowPeak six_1_peaks.narrowPeak > sorted.six1_peaks.narrowPeak


# In[8]:


head -10000 sorted.six1_peaks.narrowPeak > head10000.sorted.six1_peaks.narrowPeak


# In[9]:


cat head10000.sorted.six1_peaks.narrowPeak | cut -f 1-3 > six1_10k_peaks.bed


# In[ ]:


computeMatrix reference-point     -S six1.bw     -R six1_10k_peaks.bed     --referencePoint center     -a 2000     -b 2000     --binSize 5     -out six1_10k.tab.gz


# In[19]:


plotHeatmap     -m six1_10k.tab.gz     -out six1_10k_peaks.png     --heatmapHeight 15      --refPointLabel peak.center     --regionsLabel peaks     --plotTitle 'ChIP-seq signal'


# ### Motif Analysis

# #### 1. De novo motif discovery
# 
# Run RSAT peak-motifs at the RSAT website http://rsat.sb-roscoff.fr
# 
# Here is a tutorial for RSAT peak-motifs: http://embnet.ccg.unam.mx/rsat//htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_peak-motifs.html
# 
# Background information about position-specific scoring matrices (PSSM) is here: http://embnet.ccg.unam.mx/rsat//htmllink.cgi?title=RSAT-tutorials&file=tutorials/tut_peak-motifs.html
# 
# Before running peak-motifs, you need as input the <b>FASTA sequences of the peaks</b>. How to get these?
# One option: use fetch-sequences from UCSC in RSAT first. http://rsat.sb-roscoff.fr/fetch-sequences_form.php.
# Other options include UCSC Table Browser (first make a custom track of your peak BED file) or command-line (e.g., using BEDTools).
# 
# 

# See the master notebook for the outcome of RSAT. 

# #### Create a BED file with only the DIRECT peaks = peaks with the TP53 motif
# 
# 1253 peaks have a predicted site, you can download these by clicking on the "sites" link on the RSAT results page.

# In[1]:


cat peak-motifs_oligos_6nt_mkv4_m1_sites.tab | grep -v ";" | grep -v '#' | head


# In[2]:


cat peak-motifs_oligos_6nt_mkv4_m1_sites.tab | grep -v ";" | grep -v '#' | cut -f 1 | tr "_" "\t" | cut -f 2-4 | head


# In[5]:


cat peak-motifs_oligos_6nt_mkv4_m1_sites.tab | grep -v ";" | grep -v '#' | cut -f 1 | tr "_" "\t" | cut -f 2-4 > six1-allpeaks-with-motif-RSAT.bed


# In[6]:


cat six1-allpeaks-with-motif-RSAT.bed | wc -l


# ### PWM and track enrichment
# Next to de novo motif enrichment, which we did above with RSAT, we can also ask whether the TP53 ChIP-seq peaks we found are enriched for matches to known PSSMs (PWMs) from PSSM databases like JASPAR, TRANFAC, cis-bp, HOCOMOCO, etc.
# 
# In addition, we can ask whether our peaks overlap significantly (or "correlate") with other epigenomic data sets, called "tracks". For example, histone modification ChIP-seq tracks, or chromatin accessibility tracks. Some web-based tools that do this include i-cisTarget, Enrichr, and LOLA.
# 
# i-cisTarget will test for enrichment of motifs, using a collection of databases of PSSMs, as well as tracks.
# 
# Below, we have uploaded our ChIP-seq peak file to i-cisTarget https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/. This will place your job in a queue, which can take a while (~15 min or so) to complete.

# #### For the whole set of peaks (10.000 peaks) 

# In[12]:


head 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt


# In[13]:


cat 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt | cut -f 1 | grep -v '#' | wc -l


# In[15]:


cat 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt | cut -f 1 | grep -v '#' | grep Slc16a10
cat 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt | cut -f 1 | grep -v '#' | grep MYOG
cat 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt | cut -f 1 | grep -v '#' | grep IFGBP5

#none of the control genes is present


# In[21]:


cat 20211219-public-4.0.4-NiPCRW-mm10-all-gene.txt | cut -f 1 | grep -v '#' > Six1-10k-targets-GREAT.txt


# #### For the direct peaks (92 peaks) 

# In[ ]:


head 20211219-public-4.0.4-aCWWPy-mm10-all-gene.txt


# In[18]:


cat 20211219-public-4.0.4-aCWWPy-mm10-all-gene.txt | cut -f 1 | grep -v '#' | wc -l


# In[22]:


cat 20211219-public-4.0.4-aCWWPy-mm10-all-gene.txt | cut -f 1 | grep -v '#' > Six1-92-targets-GREAT.txt


# ### Comparison of predicted targets with functional associations to Six1

# Now, we head to the String to visualize the gene network. Let's see if some more of these genes are "related" to Six1. Using https://string-db.org/ we can obtain a network of associated genes with Six1, based on genetic interactions, physical interactions, co-expression, and literature co-mentioning (in abstracts).
# 
# 
# See the master notebook for the results and interpretation. 
