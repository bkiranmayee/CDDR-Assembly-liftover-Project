# CDDR-Assembly-liftover-Project
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
sbatch -p assemble2 faSplit_minimap2_align_kb.sh ARS-UCD1.0.14.clean.wIGCHaps.fasta.mmi /mnt/nfs/nfs2/Genomes/umd3_kary_unmask_ngap.fa
sbatch -p assemble2 faSplit_minimap2_align_kb.r2.sh umd3_kary_unmask_ngap.fa.mmi ARS-UCD1.0.14.clean.wIGCHaps.fasta

```

The jobs take around 7.5 hrs to run...so waiting
