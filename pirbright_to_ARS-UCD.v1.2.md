# Liftover of pirbright alternate haplotypes to the new assembly ARS-UCDv1.2 

  * I will be using minimap2 for genome alignment which gives a sam output file. 
  * I will then convert the sam file to psl, then chain format and finally the liftover chain

Working directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers

### First the forward direction
Old assembly (From) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc
New assembly (To) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi

```bash
sbatch -p assemble3 minimap2_align_kb.sh /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/pirbright_combined.fa /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa
 
/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl liftovers/pirbright_to_ARS-UCDv1.2.psl /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/pirbright_combined.2bit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa.2bit liftovers/pirbright_to_ARS-UCDv1.2.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort liftovers/pirbright_to_ARS-UCDv1.2.chain liftovers/pirbright_to_ARS-UCDv1.2.sorted.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet liftovers/pirbright_to_ARS-UCDv1.2.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/pirbright_combined.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa.2bit.info liftovers/pirbright_to_ARS-UCDv1.2.net /dev/null

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset liftovers/pirbright_to_ARS-UCDv1.2.net liftovers/pirbright_to_ARS-UCDv1.2.sorted.chain liftovers/pirbright_to_ARS-UCDv1.2.liftover.chain
```

### Liftover of IGC VCF files to ARS-UCDv1.2

First generate the dictionary of the reference/new assembly or else the picard LiftoverVcf will through a java error that is not very informative as I cannot understand and keep wondering what the problem was!!!

```bash
module load java
module load picard

java -jar $PICARD CreateSequenceDictionary R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa O=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa.dict

# check for the liftover filtered IGC VCF file
[kiranmayee.bakshy@assembler2 condensed_vcfs]$ java -jar $PICARD LiftoverVcf I=igc.filtered2.vcf.gz O=liftover_to_v1.2/igc.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/pirbright_to_ARS-UCDv1.2/pirbright_to_ARS-UCDv1.2.liftover.chain REJECT=liftover_to_v1.2/igc.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa WRITE_ORIGINAL_POSITION=true
```

VCF Liftover was 16.3% successful. 



