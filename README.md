# CDDR-Assembly-liftover-Project
Started on Apr 18 2018
This is a process to create a liftover chain file for converting the assembly coordinates from 
the new cattle assembly ars-ucd.v14 to the older version umd3 

The overall preocess is a general one that can be found by googling assembly liftover process...

Non-repeatmasked chain of assembly file

The following steps are already done:

Derek's comment: I need to check to see if I can get this working on my non-repeatmasked fastas for full liftover of v14 to umd3.

## working directory:
Assembler2: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc

# The query assembly fasta is first divided into chr chunks of 10000 bases for the BLAT aligner to work faster
mkdir chrchunks
perl -lane 'system("samtools faidx ARS-UCD1.0.14.clean.wIGCHaps.fasta $F[0] > chrchunks/$F[0].fa");' < ARS-UCD1.0.14.clean.wIGCHaps.fasta.fai

# Converting to 2bit assemblies
for i in chrchunks/*.fa; do name=`basename $i | cut -d'.' -f1`; echo $name; /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/faToTwoBit $i chrchunks/${name}.2bit; done

# start aligning using BLAT
module load blat
for i in chrchunks/*.fa; do sbatch -p assemble3 faSplit_blat_align.sh umd3_kary_unmask_ngap.2bit $i; done

The above script 'faSplit_blat_align.sh' generates .psl files which have to be 'chained'

I started working from this step:

# start chaining
sbatch -p assemble2 psl_merge.sh ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit umd3_kary_unmask_ngap.2bit

The script psl_merge.sh generates chains in the directory working_dir/umd3_psl/chr_sub_chunks which are then grouped per chr 
and a chain file per chr is created.

These have to be combined again using the following command:
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainMergeSort /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/*_chain/*.chain | /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSplit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/umd3_psl/combined_chain_blat stdin

We are not bothered about the unscaffolded elements because the liftover is from a new version (ars.v14) to an old version (umd3)
We mainly require liftover of the main chr (1-29 and X)

Merging all the chr chains
cat umd3_psl/combined_chain_blat/*.chain > umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.chain

Sorting the chain
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.chain umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain
 
preparing for net
mkdir umd3_psl/net
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/twoBitInfo ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info

Creating net file
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain umd3_kary_unmask_ngap.2bit.info ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_psl/net/combined.ars.v14_to_umd3.net /dev/null
Got 30 chroms in umd3_kary_unmask_ngap.2bit.info, 30 in ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info
Finishing nets
writing umd3_psl/net/combined.ars.v14_to_umd3.net
writing /dev/null

Creating the final liftover chain
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet umd3_psl/combined_chain_blat/combined.ars.v14_to_umd3.sorted.chain ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info umd3_kary_unmask_ngap.2bit.info umd3_psl/net/combined.ars.v14_to_umd3.net /dev/null
Got 2217 chroms in ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info, 30 in umd3_kary_unmask_ngap.2bit.info
hashMustFindVal: '24' not found
