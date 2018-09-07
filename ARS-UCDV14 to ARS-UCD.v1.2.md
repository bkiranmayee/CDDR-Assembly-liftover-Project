#Liftover of ARS-UCDv14 to and from the new assembly ARS-UCDv1.2 

  * I will be using minimap2 for genome alignment which gives a sam output file. 
  * I will then convert the sam file to psl, then chain format and finally the liftover chain

Working directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers

###First the forward direction
Old assembly (From) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc
New assembly (To) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi

```bash
sbatch -p assemble3 minimap2_align_kb.sh /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.mmi /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl liftovers/ARS-UCDv14_to_ARS-UCDv1.2.psl /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit.info liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net /dev/null

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain
```

####Check the liftover performance:

```bash
[kiranmayee.bakshy@assembler2 ARS-UCDv14_to_ARS-UCDv1.2]$ /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/test.bed ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain mapped unmapped
Reading liftover chains
Mapping coordinates
[kiranmayee.bakshy@assembler2 ARS-UCDv1.2_to_ARS-UCDv14]$ wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/test.bed mapped unmapped
  23912824 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/vcf2bed/test.bed
  23782454 mapped  <- 99.5%
    260740 unmapped <- actually 130370 positions i.e. 0.5%
  47956018 total
```

###Now the reverse direction
Old assembly (From) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi
New assembly (To) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc


```bash
sbatch -p assemble3 minimap2_align_kb.sh /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.mmi /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl liftovers/ARS-UCDv14_to_ARS-UCDv1.2.psl /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain 

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net /dev/null

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain
```