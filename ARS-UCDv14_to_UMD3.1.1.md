# Liftover ARS-UCDv14 to UMD3.1.1
---

*Apr 18 2018*

These are my notes to create a liftover chain file for converting the assembly coordinates from 
the new cattle assembly ars-ucd.v14 to the older version umd3 

## Table of Contents
 * [Liftover using BLAT alignment](#BLAT)
     * [Run2 using exchanged assemblies](#exchangedassemblies)
 * [Liftover using Minimap2 genome to genome alignment](#MINIMAP2)
 
<a name="BLAT"></a>
## Liftover using BLAT alignment

The overall process is a general one that can be found by googling assembly liftover process...

Non-repeatmasked chain of assembly file

The following steps were already done by Derek

> Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

```bash
# The query assembly fasta is first divided into chr chunks of 10000 bases for the BLAT aligner to work faster
mkdir chrchunks
perl -lane 'system("samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0] > chrchunks/$F[0].fa");' < ARS-UCD1.0.14.clean.wIGCHaps.fasta.fai

# Converting to 2bit assemblies
for i in chrchunks/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/faToTwoBit $i chrchunks/${name}.2bit; done

# start aligning using BLAT
module load blat
for i in chrchunks/*.fa; do sbatch -p assemble3 faSplit_blat_align.sh umd3_kary_unmask_ngap.2bit $i; done
```

The above script 'faSplit_blat_align.sh' generates .psl files which have to be 'chained'

I started working from this step:

```bash
# start chaining
sbatch -p assemble2 psl_merge.sh ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit umd3_kary_unmask_ngap.2bit

# The script psl_merge.sh generates chains in the directory working_dir/umd3_psl/chr_sub_chunks which are then grouped per chr 
# and a chain file per chr is created.

# These have to be combined again using the following command:
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/*_chain/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat stdin
```

We are not bothered about the unscaffolded elements because the liftover is from a new version (ars.v14) to an old version (umd3)
We mainly require liftover of the main chr (1-29 and X)

```bash
# Merging all the chr chains
cat umd3_psl/combined_chain_blat/*.chain > umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.chain

# Sorting the chain
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.chain umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain
 
# preparing for net
mkdir umd3_psl/net
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/twoBitInfo ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info

# Creating net file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_kary_unmask_ngap.2bit.info umd3_psl/net/combined.ars.v14_to_umd3.net /dev/null
Got 2217 chroms in ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info, 30 in umd3_kary_unmask_ngap.2bit.info
hashMustFindVal: '24' not found

```
The error suggests that the query and ref fastas should be reversed. 
Quickly checking the fasta headers of the files by using the samtools “faidx” command.
From the fai files, the chromosome “24” is present in both the files 

However, after reversing the query and ref for the chainNet command, there is no error and the net file is ready
```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain umd3_kary_unmask_ngap.2bit.info ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_psl/net/combined.ars.v14_to_umd3.net /dev/null
Got 30 chroms in umd3_kary_unmask_ngap.2bit.info, 2217 in ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info
Finishing nets
writing umd3_psl/net/combined.ars.v14_to_umd3.net
writing /dev/null

#Now the final step:
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset umd3_psl/net/combined.ars.v14_to_umd3.net umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain umd3_psl/net/combined.ars.v14_to_umd3.liftover.chain
```

Finally check the liftover process by lifting the bed coordinates from the variant calling project

> Assembler2: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs

Make bed coordinate files for all the variants in these vcf files
```bash
mkdir vcf2bed
module load bcftools bedtools
for i in `seq 1 29` X; do bcftools query -f 'chr%CHROM\t%POS0\t%END\n'  ${i}.vcf.gz -o vcf2bed/${i}.ars.v14.bed; done

# merge all the bed files
cat vcf2bed/*.ars.v14.bed > vcf2bed/var_ars_ucd.v14.bed

# sort bed file
sortBed -chrThenSizeA -i vcf2bed/var_ars_ucd.v14.bed > vcf2bed/var_ars_ucd.v14.sorted.bed

# Now do the liftover 
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver vcf2bed/var_ars_ucd.v14.sorted.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/combined.ars.v14_to_umd3.liftover.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var_umd3.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var.umd3.bed.unmapped
 
 # And check
 wc -l vcf2bed/var_ars_ucd.v14.sorted.bed  /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var_umd3.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var.umd3.bed.unmapped
 23912824 vcf2bed/var_ars_ucd.v14.sorted.bed
  1040884 /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var_umd3.bed    <- 4.35% are mapped to umd3 which is the opposite of what we should get
 45743880 /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/net/var.umd3.bed.unmapped   <- 96% unmapped
```

It needs further interrogation

Inspect .psl files and .chain files

Check the steps and the files generated.


<a name="exchangedassemblies"></a>
#### Run2 using exchanged assemblies

Trying to run the psl merge step with assemblies exchanged for ref and query:

Script was modified to generate the new chain files with R2 suffix that indicates run2 (psl_merge_kb.sh)

Working directory : /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl

```bash
# start chaining
sbatch -p assemble2 psl_merge_kb.sh umd3_kary_unmask_ngap.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit

# The script worked partially and failed to generate merged chains, the slurm-815776.out says 
Permission denied
Couldn't make directory /X_chain.R2
But the directories were actually created which had no files in it

# Merging these new chain files in a separate step
for i in `seq 1 29` X; do /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/chr_sub_chunks/${i}_*.R2.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/${i}_R2.chain stdin; done

# Combining the chains using the following command:
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/*_R2.chain/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2 stdin

# Merging all the chr chains
cat umd3_psl/combined_chain_blat.R2/*.chain > umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.chain

# Sorting the chain
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.chain umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.sorted.chain
 
# preparing for net
mkdir umd3_psl/net
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/twoBitInfo ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info

# Creating net file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.sorted.chain umd3_kary_unmask_ngap.2bit.info ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.net /dev/null
Got 30 chroms in umd3_kary_unmask_ngap.2bit.info, 30 in ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info
Finishing nets
writing umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.net
writing /dev/null

#Now the final step:
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.net umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.sorted.chain umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.liftover.chain
Processing chr1
Processing chr2
Processing chr3
Processing chr4
Processing chr5
Processing chr6
Processing chr7
Processing chr8
Processing chr9
Processing chr10
Processing chr11
Processing chr12
Processing chr13
Processing chr14
Processing chr15
Processing chr16
Processing chr17
Processing chr18
Processing chr19
Processing chr20
Processing chr21
Processing chr22
Processing chr23
Processing chr24
Processing chr25
Processing chr26
Processing chr27
Processing chr28
Processing chr29
Processing chrX
```

OK, now the new liftover file is ready, I would like to check it.

```bash
grep 'chain' /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.liftover.chain | perl -lane 'print $F[2];' | perl perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -m

The above command didn't work, some PATH issues or something else 
nameListVennCount.pl script from the same folder works well but tabFileColumnCounter.pl does not.
says: Can't locate Mouse.pm

# doing the liftover
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/combined.ars.v14_to_umd3.liftover.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var_umd3.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var.umd3.bed.unmapped

# Check
 wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed  /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var_umd3.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var.umd3.bed.unmapped
  23912824 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed
  22041746 /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var_umd3.bed
   3742156 /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat.R2/var.umd3.bed.unmapped <- acutally 1871078 unmapped coordinates
  49696726 total
  mapped 92.2 %
  unmapped 7.8%
  ```
  Looks good!
  
There are only 3 comments in the unmapped file: 
 * Deleted in new, 
 * Partially deleted in new,
 * Split in new.

Checking no. of mapped and unmapped positions per chr
```bash
for i in `seq 1 29` X; do grep -P -c "^${i}\t" var_umd3.bed; done
for i in `seq 1 29` X; do grep -P -c "^chr${i}\t" var.umd3.bed.unmapped; done

```
| Chr   | unmapped | mapped   |
|-------|----------|----------|
| 1     | 50452    | 1446638  |
| 2     | 26812    | 1172818  |
| 3     | 38022    | 1020819  |
| 4     | 33378    | 1107433  |
| 5     | 35349    | 1061161  |
| 6     | 89493    | 1009999  |
| 7     | 49568    | 909069   |
| 8     | 40039    | 940160   |
| 9     | 36039    | 935833   |
| 10    | 69117    | 877597   |
| 11    | 23048    | 877464   |
| 12    | 141361   | 857306   |
| 13    | 23416    | 623273   |
| 14    | 28619    | 681012   |
| 15    | 42354    | 814693   |
| 16    | 40104    | 638619   |
| 17    | 34816    | 672204   |
| 18    | 48984    | 546515   |
| 19    | 20948    | 490303   |
| 20    | 20973    | 653529   |
| 21    | 28220    | 583328   |
| 22    | 13975    | 469894   |
| 23    | 65575    | 573192   |
| 24    | 14962    | 543848   |
| 25    | 13234    | 360100   |
| 26    | 13538    | 450934   |
| 27    | 31771    | 429846   |
| 28    | 16412    | 433489   |
| 29    | 31208    | 498538   |
| X     | 749291   | 362132   |
| total | 1871078  | 22041746 |


<a name="MINIMAP2"></a>
## Liftover using MINIMAP2 genome to genome alignment

* Whole genome sequence alignment using minimap2 is much faster than that using BLAT alignment. 
* The process is simple first use minimap2 for aligning the genomes fasta files which outputs sam format.
* The sam format can be converted into psl using a python converter program. This step is important to compare the liftover methods.

    The minimap2 aligment uses minimizer files (.mmi) which is a kind of index for the reference fasta file.
    
    So first generate .mmi for the reference and use it for the batch job.
    
```bash
job1 sbatch -p assemble2 faSplit_minimap2_align_kb.sh ARS-UCD1.0.14.clean.wIGCHaps.fasta.mmi /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa
job2 sbatch -p assemble2 faSplit_minimap2_align_kb.r2.sh umd3_kary_unmask_ngap.fa.mmi ARS-UCD1.0.14.clean.wIGCHaps.fasta

```

The jobs take around 7.5 hrs to run...so waiting
Actually job2 was completed in 20 min!!!

Now making chain from psl

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl minimap_liftover/umd3_ars.v14_mmap.r2.psl /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa ARS-UCD1.0.14.clean.wIGCHaps.fasta minimap_liftover/umd3_ars.v14_mmap.r2.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl minimap_liftover/umd3_ars.v14_mmap.r2.psl ARS-UCD1.0.14.clean.wIGCHaps.fasta /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa minimap_liftover/umd3_ars.v14_mmap.r2.chain

showing an error 
invalid unsigned integer: "-189"
Line 39813 of the psl file has the first field (Number of matching bases that aren't repeats) -189 
I changed it to 189 and lets see if it works now?
Next the error is different:
it says Can't open query fasta because it is not a directory
Upon looking at the axtChain script looks like I should be giving the 2bit files for query and target instead of fasta

So here it is:

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl minimap_liftover/umd3_ars.v14_mmap.r2.psl umd3_kary_unmask_ngap.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit minimap_liftover/umd3_ars.v14_mmap.r2.chain

Yes! Now it is working...

Sorting chain...
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/minimap_liftover/umd3_ars.v14_mmap.r2.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/minimap_liftover/umd3_ars.v14_mmap.r2.sorted.chain

Generating net
 /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet minimap_liftover/umd3_ars.v14_mmap.r2.sorted.chain umd3_kary_unmask_ngap.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info minimap_liftover/umd3_ars.v14_mmap.r2.net /dev/null

Got 30 chroms in umd3_kary_unmask_ngap.2bit.info, 2217 in /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info
Finishing nets
writing minimap_liftover/umd3_ars.v14_mmap.r2.net
writing /dev/null

Finally the liftover chain file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset minimap_liftover/umd3_ars.v14_mmap.r2.net  minimap_liftover/umd3_ars.v14_mmap.r2.sorted.chain minimap_liftover/umd3_ars.v14_mmap.r2.liftover.chain
Processing chr1
Processing chr2
Processing chr3
Processing chr4
Processing chr5
Processing chr6
Processing chr7
Processing chr8
Processing chr9
Processing chr10
Processing chr11
Processing chr12
Processing chr13
Processing chr14
Processing chr15
Processing chr16
Processing chr17
Processing chr18
Processing chr19
Processing chr20
Processing chr21
Processing chr22
Processing chr23
Processing chr24
Processing chr25
Processing chr26
Processing chr27
Processing chr28
Processing chr29
Processing chrX

Suuccessfully completed...
```
Now check the liftover...

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed minimap_liftover/umd3_ars.v14_mmap.r2.liftover.chain minimap_liftover/var_umd3.bed minimap_liftover/var.umd3.bed.unmapped

 wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed minimap_liftover/var_umd3.bed minimap_liftover/var.umd3.bed.unmapped
  23912824 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed
  23145319 minimap_liftover/var_umd3.bed <- 96.8%
   1535010 minimap_liftover/var.umd3.bed.unmapped  <- actually 767505 unmapped positions i.e. 3.2%
  48593153 total
   ```
  Looks great!
  
There are only 3 comments in the unmapped file: 
 * Deleted in new, 
 * Partially deleted in new,
 * Split in new.

Checking no. of mapped and unmapped positions per chr
```bash
for i in `seq 1 29` X; do grep -P -c "^${i}\t" minimap_liftover/var_umd3.bed; done
for i in `seq 1 29` X; do grep -P -c "^chr${i}\t" minimap_liftover/var.umd3.bed.unmapped; done
```
| Chr         | unmapped | mapped   |
|-------------|----------|----------|
| 1           | 26559    | 1469518  |
| 2           | 12274    | 1185709  |
| 3           | 16370    | 1041343  |
| 4           | 17698    | 1121037  |
| 5           | 16615    | 1076721  |
| 6           | 75886    | 1021407  |
| 7           | 28140    | 927566   |
| 8           | 19265    | 956678   |
| 9           | 18954    | 951118   |
| 10          | 30467    | 911119   |
| 11          | 10982    | 888579   |
| 12          | 37149    | 905238   |
| 13          | 12790    | 630078   |
| 14          | 16447    | 689507   |
| 15          | 21041    | 834161   |
| 16          | 15367    | 655040   |
| 17          | 23834    | 681125   |
| 18          | 19495    | 578807   |
| 19          | 9369     | 500211   |
| 20          | 9838     | 663993   |
| 21          | 16053    | 593359   |
| 22          | 10456    | 473403   |
| 23          | 39496    | 598829   |
| 24          | 6032     | 547886   |
| 25          | 7328     | 365287   |
| 26          | 8428     | 455645   |
| 27          | 15145    | 442221   |
| 28          | 6151     | 442868   |
| 29          | 10718    | 514160   |
| X           | 209158   | 882014   |
| total       | 767505   | 23004627 |
| Leftovers   |          | 140692   |
| Grand total |          | 23145319 |


I have to complete my run 1...
So,
```bash
mkdir minimap_liftover/r1
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl minimap_liftover/r1/umd3_ars.v14_mmap.psl ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit umd3_kary_unmask_ngap.2bit minimap_liftover/r1/umd3_ars.v14_mmap.chain
invalid unsigned integer: "-88"
invalid unsigned integer: "-276" line no. 3090
invalid unsigned integer: "-11985" line no. 6798
invalid unsigned integer: "-2079"   line no. 14385
I edited these to positive values in the file
Yes! Now it is working...

Sorting chain...
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/minimap_liftover/r1/umd3_ars.v14_mmap.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/minimap_liftover/r1/umd3_ars.v14_mmap.sorted.chain

Generating net
 /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet minimap_liftover/r1/umd3_ars.v14_mmap.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_kary_unmask_ngap.2bit.info minimap_liftover/r1/umd3_ars.v14_mmap.net /dev/null
 Got 2217 chroms in /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info, 30 in umd3_kary_unmask_ngap.2bit.info
Finishing nets
writing minimap_liftover/r1/umd3_ars.v14_mmap.net
writing /dev/null


Finally the liftover chain file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset minimap_liftover/r1/umd3_ars.v14_mmap.net  minimap_liftover/r1/umd3_ars.v14_mmap.sorted.chain minimap_liftover/r1/umd3_ars.v14_mmap.liftover.chain
Processing 1
Processing 2
Processing 3
Processing 4
Processing 5
Processing 6
Processing 7
Processing 8
Processing 9
Processing 10
Processing 11
Processing 12
Processing 13
Processing 14
Processing 15
Processing 16
Processing 17
Processing 18
Processing 19
Processing 20
Processing 21
Processing 22
Processing 23
Processing 24
Processing 25
Processing 26
Processing 27
Processing 28
Processing 29
Processing X
Processing contig_X_unplaced
Processing Leftover_ScbfJmS_100
Processing Leftover_ScbfJmS_1014
Processing Leftover_ScbfJmS_1015
Processing Leftover_ScbfJmS_1020
Processing Leftover_ScbfJmS_1022
Processing Leftover_ScbfJmS_1023
Processing Leftover_ScbfJmS_102
Processing Leftover_ScbfJmS_1044
Processing Leftover_ScbfJmS_1045
Processing Leftover_ScbfJmS_1052
Processing Leftover_ScbfJmS_1072
Processing Leftover_ScbfJmS_1074
Processing Leftover_ScbfJmS_1083
Processing Leftover_ScbfJmS_1085
Processing Leftover_ScbfJmS_109
Processing Leftover_ScbfJmS_1095
Processing Leftover_ScbfJmS_1096
Processing Leftover_ScbfJmS_1097
Processing Leftover_ScbfJmS_1107_221
Processing Leftover_ScbfJmS_111
Processing Leftover_ScbfJmS_1120
Processing Leftover_ScbfJmS_1136
Processing Leftover_ScbfJmS_1151
Processing Leftover_ScbfJmS_115
Processing Leftover_ScbfJmS_1153
Processing Leftover_ScbfJmS_116
Processing Leftover_ScbfJmS_1170_874
Processing Leftover_ScbfJmS_1177
Processing Leftover_ScbfJmS_1180
Processing Leftover_ScbfJmS_1187
Processing Leftover_ScbfJmS_120
Processing Leftover_ScbfJmS_1216_1546_2082
Processing Leftover_ScbfJmS_1218
Processing Leftover_ScbfJmS_1220
Processing Leftover_ScbfJmS_1221
Processing Leftover_ScbfJmS_1223
Processing Leftover_ScbfJmS_1233
Processing Leftover_ScbfJmS_1236
Processing Leftover_ScbfJmS_1239
Processing Leftover_ScbfJmS_1241
Processing Leftover_ScbfJmS_1248
Processing Leftover_ScbfJmS_1260
Processing Leftover_ScbfJmS_126
Processing Leftover_ScbfJmS_1270
Processing Leftover_ScbfJmS_1274
Processing Leftover_ScbfJmS_1280
Processing Leftover_ScbfJmS_1294
Processing Leftover_ScbfJmS_1296
Processing Leftover_ScbfJmS_1303
Processing Leftover_ScbfJmS_1307
Processing Leftover_ScbfJmS_1309
Processing Leftover_ScbfJmS_1312
Processing Leftover_ScbfJmS_1316
Processing Leftover_ScbfJmS_1325
Processing Leftover_ScbfJmS_1329
Processing Leftover_ScbfJmS_1330
Processing Leftover_ScbfJmS_1336
Processing Leftover_ScbfJmS_1343
Processing Leftover_ScbfJmS_1346
Processing Leftover_ScbfJmS_1349
Processing Leftover_ScbfJmS_1350
Processing Leftover_ScbfJmS_1357
Processing Leftover_ScbfJmS_1361
Processing Leftover_ScbfJmS_1368
Processing Leftover_ScbfJmS_1375
Processing Leftover_ScbfJmS_1377
Processing Leftover_ScbfJmS_1382
Processing Leftover_ScbfJmS_1386
Processing Leftover_ScbfJmS_1388
Processing Leftover_ScbfJmS_1391
Processing Leftover_ScbfJmS_1394
Processing Leftover_ScbfJmS_1396
Processing Leftover_ScbfJmS_1398
Processing Leftover_ScbfJmS_1402
Processing Leftover_ScbfJmS_1404
Processing Leftover_ScbfJmS_1409
Processing Leftover_ScbfJmS_1411
Processing Leftover_ScbfJmS_1413
Processing Leftover_ScbfJmS_1415
Processing Leftover_ScbfJmS_1416
Processing Leftover_ScbfJmS_1417
Processing Leftover_ScbfJmS_1420
Processing Leftover_ScbfJmS_1424
Processing Leftover_ScbfJmS_1427
Processing Leftover_ScbfJmS_1431
Processing Leftover_ScbfJmS_1432
Processing Leftover_ScbfJmS_143
Processing Leftover_ScbfJmS_1433
Processing Leftover_ScbfJmS_1437
Processing Leftover_ScbfJmS_14
Processing Leftover_ScbfJmS_1441
Processing Leftover_ScbfJmS_1450
Processing Leftover_ScbfJmS_1451_926
Processing Leftover_ScbfJmS_1456
Processing Leftover_ScbfJmS_1457
Processing Leftover_ScbfJmS_146
Processing Leftover_ScbfJmS_1465
Processing Leftover_ScbfJmS_1466
Processing Leftover_ScbfJmS_1470
Processing Leftover_ScbfJmS_1488
Processing Leftover_ScbfJmS_149
Processing Leftover_ScbfJmS_1493
Processing Leftover_ScbfJmS_1505_733
Processing Leftover_ScbfJmS_1506
Processing Leftover_ScbfJmS_151
Processing Leftover_ScbfJmS_1511_869
Processing Leftover_ScbfJmS_1513_2139
Processing Leftover_ScbfJmS_1517
Processing Leftover_ScbfJmS_1524
Processing Leftover_ScbfJmS_1554
Processing Leftover_ScbfJmS_1565
Processing Leftover_ScbfJmS_1569
Processing Leftover_ScbfJmS_1574
Processing Leftover_ScbfJmS_1578
Processing Leftover_ScbfJmS_1580
Processing Leftover_ScbfJmS_1586
Processing Leftover_ScbfJmS_1588
Processing Leftover_ScbfJmS_159
Processing Leftover_ScbfJmS_1593
Processing Leftover_ScbfJmS_1616
Processing Leftover_ScbfJmS_1618
Processing Leftover_ScbfJmS_162
Processing Leftover_ScbfJmS_1635
Processing Leftover_ScbfJmS_1639
Processing Leftover_ScbfJmS_163
Processing Leftover_ScbfJmS_1641
Processing Leftover_ScbfJmS_1648
Processing Leftover_ScbfJmS_166
Processing Leftover_ScbfJmS_1665
Processing Leftover_ScbfJmS_1668
Processing Leftover_ScbfJmS_1671
Processing Leftover_ScbfJmS_1678
Processing Leftover_ScbfJmS_1689
Processing Leftover_ScbfJmS_1690
Processing Leftover_ScbfJmS_169
Processing Leftover_ScbfJmS_1696
Processing Leftover_ScbfJmS_1707
Processing Leftover_ScbfJmS_1709
Processing Leftover_ScbfJmS_1712
Processing Leftover_ScbfJmS_1713
Processing Leftover_ScbfJmS_1716
Processing Leftover_ScbfJmS_1728
Processing Leftover_ScbfJmS_17
Processing Leftover_ScbfJmS_1737
Processing Leftover_ScbfJmS_1748
Processing Leftover_ScbfJmS_1751
Processing Leftover_ScbfJmS_1758
Processing Leftover_ScbfJmS_1760
Processing Leftover_ScbfJmS_1767
Processing Leftover_ScbfJmS_1771
Processing Leftover_ScbfJmS_1772
Processing Leftover_ScbfJmS_1778
Processing Leftover_ScbfJmS_1784
Processing Leftover_ScbfJmS_1791
Processing Leftover_ScbfJmS_1795
Processing Leftover_ScbfJmS_1803
Processing Leftover_ScbfJmS_1814
Processing Leftover_ScbfJmS_1815
Processing Leftover_ScbfJmS_1821
Processing Leftover_ScbfJmS_1828
Processing Leftover_ScbfJmS_18
Processing Leftover_ScbfJmS_1849
Processing Leftover_ScbfJmS_1853
Processing Leftover_ScbfJmS_1856
Processing Leftover_ScbfJmS_1863
Processing Leftover_ScbfJmS_1869
Processing Leftover_ScbfJmS_1877
Processing Leftover_ScbfJmS_1880
Processing Leftover_ScbfJmS_1888
Processing Leftover_ScbfJmS_1893
Processing Leftover_ScbfJmS_1896
Processing Leftover_ScbfJmS_1899
Processing Leftover_ScbfJmS_1901
Processing Leftover_ScbfJmS_1905
Processing Leftover_ScbfJmS_1929
Processing Leftover_ScbfJmS_1941
Processing Leftover_ScbfJmS_1953
Processing Leftover_ScbfJmS_1959
Processing Leftover_ScbfJmS_1962
Processing Leftover_ScbfJmS_1967
Processing Leftover_ScbfJmS_1978
Processing Leftover_ScbfJmS_1989_199
Processing Leftover_ScbfJmS_1995
Processing Leftover_ScbfJmS_1996
Processing Leftover_ScbfJmS_2002
Processing Leftover_ScbfJmS_2005
Processing Leftover_ScbfJmS_2006
Processing Leftover_ScbfJmS_2007
Processing Leftover_ScbfJmS_2012
Processing Leftover_ScbfJmS_2017
Processing Leftover_ScbfJmS_2018
Processing Leftover_ScbfJmS_2021
Processing Leftover_ScbfJmS_202
Processing Leftover_ScbfJmS_2032
Processing Leftover_ScbfJmS_2048
Processing Leftover_ScbfJmS_2049
Processing Leftover_ScbfJmS_2053
Processing Leftover_ScbfJmS_2057
Processing Leftover_ScbfJmS_2063
Processing Leftover_ScbfJmS_2066
Processing Leftover_ScbfJmS_2072
Processing Leftover_ScbfJmS_2074_122
Processing Leftover_ScbfJmS_2090
Processing Leftover_ScbfJmS_2093
Processing Leftover_ScbfJmS_2095
Processing Leftover_ScbfJmS_2101
Processing Leftover_ScbfJmS_210
Processing Leftover_ScbfJmS_2110
Processing Leftover_ScbfJmS_2111
Processing Leftover_ScbfJmS_2116
Processing Leftover_ScbfJmS_212
Processing Leftover_ScbfJmS_2123
Processing Leftover_ScbfJmS_2125
Processing Leftover_ScbfJmS_2127
Processing Leftover_ScbfJmS_2134
Processing Leftover_ScbfJmS_2140
Processing Leftover_ScbfJmS_2142
Processing Leftover_ScbfJmS_2154
Processing Leftover_ScbfJmS_2157
Processing Leftover_ScbfJmS_2158
Processing Leftover_ScbfJmS_2159
Processing Leftover_ScbfJmS_2162
Processing Leftover_ScbfJmS_216
Processing Leftover_ScbfJmS_2166
Processing Leftover_ScbfJmS_2167
Processing Leftover_ScbfJmS_2168
Processing Leftover_ScbfJmS_2169
Processing Leftover_ScbfJmS_2175
Processing Leftover_ScbfJmS_2177
Processing Leftover_ScbfJmS_2182_1813_891_820
Processing Leftover_ScbfJmS_2183
Processing Leftover_ScbfJmS_2184
Processing Leftover_ScbfJmS_2185
Processing Leftover_ScbfJmS_2189
Processing Leftover_ScbfJmS_2217
Processing Leftover_ScbfJmS_2219
Processing Leftover_ScbfJmS_2222
Processing Leftover_ScbfJmS_2231
Processing Leftover_ScbfJmS_2234_990
Processing Leftover_ScbfJmS_2239
Processing Leftover_ScbfJmS_2246
Processing Leftover_ScbfJmS_2247_2147_1698
Processing Leftover_ScbfJmS_2249
Processing Leftover_ScbfJmS_225
Processing Leftover_ScbfJmS_2255
Processing Leftover_ScbfJmS_2264
Processing Leftover_ScbfJmS_2266
Processing Leftover_ScbfJmS_2281
Processing Leftover_ScbfJmS_2295
Processing Leftover_ScbfJmS_2296
Processing Leftover_ScbfJmS_2297
Processing Leftover_ScbfJmS_2299_ScbfJmS_836
Processing Leftover_ScbfJmS_2301
Processing Leftover_ScbfJmS_2306
Processing Leftover_ScbfJmS_2311
Processing Leftover_ScbfJmS_2317
Processing Leftover_ScbfJmS_2319
Processing Leftover_ScbfJmS_2322
Processing Leftover_ScbfJmS_2332
Processing Leftover_ScbfJmS_2335
Processing Leftover_ScbfJmS_2337
Processing Leftover_ScbfJmS_2341
Processing Leftover_ScbfJmS_234
Processing Leftover_ScbfJmS_2345
Processing Leftover_ScbfJmS_235
Processing Leftover_ScbfJmS_2357
Processing Leftover_ScbfJmS_2362
Processing Leftover_ScbfJmS_2364
Processing Leftover_ScbfJmS_2373
Processing Leftover_ScbfJmS_2378
Processing Leftover_ScbfJmS_2381
Processing Leftover_ScbfJmS_2390
Processing Leftover_ScbfJmS_2396
Processing Leftover_ScbfJmS_24
Processing Leftover_ScbfJmS_2420
Processing Leftover_ScbfJmS_2431
Processing Leftover_ScbfJmS_2439
Processing Leftover_ScbfJmS_2460
Processing Leftover_ScbfJmS_2461
Processing Leftover_ScbfJmS_2468_375
Processing Leftover_ScbfJmS_2469
Processing Leftover_ScbfJmS_2470
Processing Leftover_ScbfJmS_247
Processing Leftover_ScbfJmS_2474
Processing Leftover_ScbfJmS_261_1982
Processing Leftover_ScbfJmS_264
Processing Leftover_ScbfJmS_275
Processing Leftover_ScbfJmS_28
Processing Leftover_ScbfJmS_285
Processing Leftover_ScbfJmS_295
Processing Leftover_ScbfJmS_300
Processing Leftover_ScbfJmS_307
Processing Leftover_ScbfJmS_31
Processing Leftover_ScbfJmS_315_1043
Processing Leftover_ScbfJmS_317
Processing Leftover_ScbfJmS_3
Processing Leftover_ScbfJmS_323
Processing Leftover_ScbfJmS_324_987_2391
Processing Leftover_ScbfJmS_328
Processing Leftover_ScbfJmS_32
Processing Leftover_ScbfJmS_331
Processing Leftover_ScbfJmS_337
Processing Leftover_ScbfJmS_338
Processing Leftover_ScbfJmS_342
Processing Leftover_ScbfJmS_344
Processing Leftover_ScbfJmS_357
Processing Leftover_ScbfJmS_382
Processing Leftover_ScbfJmS_40
Processing Leftover_ScbfJmS_403
Processing Leftover_ScbfJmS_406
Processing Leftover_ScbfJmS_407
Processing Leftover_ScbfJmS_413
Processing Leftover_ScbfJmS_417
Processing Leftover_ScbfJmS_419
Processing Leftover_ScbfJmS_4
Processing Leftover_ScbfJmS_434
Processing Leftover_ScbfJmS_44
Processing Leftover_ScbfJmS_45
Processing Leftover_ScbfJmS_459
Processing Leftover_ScbfJmS_466
Processing Leftover_ScbfJmS_468
Processing Leftover_ScbfJmS_473
Processing Leftover_ScbfJmS_479
Processing Leftover_ScbfJmS_483
Processing Leftover_ScbfJmS_501
Processing Leftover_ScbfJmS_504
Processing Leftover_ScbfJmS_50
Processing Leftover_ScbfJmS_507
Processing Leftover_ScbfJmS_510
Processing Leftover_ScbfJmS_512_38
Processing Leftover_ScbfJmS_514
Processing Leftover_ScbfJmS_519
Processing Leftover_ScbfJmS_52
Processing Leftover_ScbfJmS_523
Processing Leftover_ScbfJmS_524
Processing Leftover_ScbfJmS_527
Processing Leftover_ScbfJmS_540
Processing Leftover_ScbfJmS_542
Processing Leftover_ScbfJmS_548
Processing Leftover_ScbfJmS_550_1916
Processing Leftover_ScbfJmS_552
Processing Leftover_ScbfJmS_556
Processing Leftover_ScbfJmS_563
Processing Leftover_ScbfJmS_567
Processing Leftover_ScbfJmS_571
Processing Leftover_ScbfJmS_591
Processing Leftover_ScbfJmS_592
Processing Leftover_ScbfJmS_594
Processing Leftover_ScbfJmS_596_672
Processing Leftover_ScbfJmS_600
Processing Leftover_ScbfJmS_60
Processing Leftover_ScbfJmS_604
Processing Leftover_ScbfJmS_605
Processing Leftover_ScbfJmS_608
Processing Leftover_ScbfJmS_610
Processing Leftover_ScbfJmS_613
Processing Leftover_ScbfJmS_626
Processing Leftover_ScbfJmS_636
Processing Leftover_ScbfJmS_637
Processing Leftover_ScbfJmS_640
Processing Leftover_ScbfJmS_646
Processing Leftover_ScbfJmS_647
Processing Leftover_ScbfJmS_659
Processing Leftover_ScbfJmS_668_1988
Processing Leftover_ScbfJmS_698
Processing Leftover_ScbfJmS_699
Processing Leftover_ScbfJmS_700
Processing Leftover_ScbfJmS_703
Processing Leftover_ScbfJmS_709
Processing Leftover_ScbfJmS_712
Processing Leftover_ScbfJmS_713
Processing Leftover_ScbfJmS_715
Processing Leftover_ScbfJmS_717
Processing Leftover_ScbfJmS_719
Processing Leftover_ScbfJmS_734
Processing Leftover_ScbfJmS_739
Processing Leftover_ScbfJmS_740
Processing Leftover_ScbfJmS_751_609
Processing Leftover_ScbfJmS_752
Processing Leftover_ScbfJmS_75
Processing Leftover_ScbfJmS_758
Processing Leftover_ScbfJmS_759
Processing Leftover_ScbfJmS_761
Processing Leftover_ScbfJmS_76
Processing Leftover_ScbfJmS_765
Processing Leftover_ScbfJmS_785
Processing Leftover_ScbfJmS_78
Processing Leftover_ScbfJmS_792
Processing Leftover_ScbfJmS_800
Processing Leftover_ScbfJmS_807
Processing Leftover_ScbfJmS_814
Processing Leftover_ScbfJmS_818
Processing Leftover_ScbfJmS_822
Processing Leftover_ScbfJmS_826
Processing Leftover_ScbfJmS_828
Processing Leftover_ScbfJmS_829
Processing Leftover_ScbfJmS_830
Processing Leftover_ScbfJmS_832
Processing Leftover_ScbfJmS_834
Processing Leftover_ScbfJmS_842
Processing Leftover_ScbfJmS_852
Processing Leftover_ScbfJmS_855
Processing Leftover_ScbfJmS_862
Processing Leftover_ScbfJmS_870
Processing Leftover_ScbfJmS_882
Processing Leftover_ScbfJmS_884_780
Processing Leftover_ScbfJmS_886
Processing Leftover_ScbfJmS_889
Processing Leftover_ScbfJmS_896
Processing Leftover_ScbfJmS_89
Processing Leftover_ScbfJmS_90
Processing Leftover_ScbfJmS_905
Processing Leftover_ScbfJmS_908_1359
Processing Leftover_ScbfJmS_91
Processing Leftover_ScbfJmS_933
Processing Leftover_ScbfJmS_935
Processing Leftover_ScbfJmS_939
Processing Leftover_ScbfJmS_940
Processing Leftover_ScbfJmS_942
Processing Leftover_ScbfJmS_946
Processing Leftover_ScbfJmS_949
Processing Leftover_ScbfJmS_954
Processing Leftover_ScbfJmS_965_547
Processing Leftover_ScbfJmS_969
Processing Leftover_ScbfJmS_972
Processing Leftover_ScbfJmS_974
Processing Leftover_ScbfJmS_975
Processing Leftover_ScbfJmS_977
Processing Leftover_ScbfJmS_998
Processing ScbfJmS_2447
Processing ScbfJmS_511
Processing ScbfJmS_708
Processing ScbfJmS_848_142
Processing ScbfJmS_892_Leftover_ScbfJmS_927
Processing Super-Scaffold_1723_ScbfJmS_2085
Processing CH240_391K10_KIR
Processing Domino_MHCclassI_gene2-5hap1_MHC
Processing CH240_370M3_LILR_LRC
Processing TPI4222_A14_MHCclassI_MHC
Processing HF_LRC_hap1_KIR_LRC
Processing LIB14427_MHC
Processing LIB14413_LRC

Suuccessfully completed...
```
Now check the liftover...

```bash
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed minimap_liftover/r1/umd3_ars.v14_mmap.liftover.chain minimap_liftover/r1/var_umd3.bed minimap_liftover/r1/var.umd3.bed.unmapped
Reading liftover chains
Mapping coordinates

 wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed minimap_liftover/r1/var_umd3.bed minimap_liftover/r1/var.umd3.bed.unmapped
    23912824 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed
         0 minimap_liftover/r1/var_umd3.bed
  47825648 minimap_liftover/r1/var.umd3.bed.unmapped
  71738472 total
   ```
  Liftover failed, there is something wrong and I guess it is with the chr notation in the bed file
 So, first change the notation from chr1 to 1 and mapp again
 
 ```bash
awk '{gsub(/^chr/,""); print}' var_ars_ucd.v14.sorted.bed > test.bed

# Lifting over...
[kiranmayee.bakshy@assembler2 ars_ucd_114_igc]$ /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/test.bed minimap_liftover/r1/umd3_ars.v14_mmap.liftover.chain minimap_liftover/r1/test_var_umd3.bed minimap_liftover/r1/test_var.umd3.bed.unmapped
Reading liftover chains
Mapping coordinates

# Checking liftover
[kiranmayee.bakshy@assembler2 ars_ucd_114_igc]$ wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed minimap_liftover/r1/test_var_umd3.bed minimap_liftover/r1/test_var.umd3.bed.unmapped
  23912824 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/var_ars_ucd.v14.sorted.bed
  23203923 minimap_liftover/r1/test_var_umd3.bed
   1417802 minimap_liftover/r1/test_var.umd3.bed.unmapped <- actually 708901 unmapped positions
  48534549 total
 
 Mapped 97%
 unmapped 3%
 
```
Looks wonderful and correct!!!

```bash
[kiranmayee.bakshy@assembler2 ars_ucd_114_igc]$ for i in `seq 1 29` X; do grep -P -c "^chr${i}\t" minimap_liftover/r1/test_var_umd3.bed; done
[kiranmayee.bakshy@assembler2 ars_ucd_114_igc]$ for i in `seq 1 29` X; do grep -P -c "^${i}\t" minimap_liftover/r1/test_var.umd3.bed.unmapped; done




  




  
