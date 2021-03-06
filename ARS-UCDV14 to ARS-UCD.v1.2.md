# Liftover of ARS-UCDv14 to and from the new assembly ARS-UCDv1.2 

  * I will be using minimap2 for genome alignment which gives a sam output file. 
  * I will then convert the sam file to psl, then chain format and finally the liftover chain

Working directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers

### First the forward direction
Old assembly (From) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc
New assembly (To) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi

```bash
sbatch -p assemble3 minimap2_align_kb.sh /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.mmi /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl liftovers/ARS-UCDv14_to_ARS-UCDv1.2.psl /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit.info liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net /dev/null

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain
```

#### Check the liftover performance:

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

### Now the reverse direction
Old assembly (From) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi
New assembly (To) directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc


```bash
sbatch -p assemble3 minimap2_align_kb.sh /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.mmi /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/axtChain -linearGap=medium -psl liftovers/ARS-UCDv14_to_ARS-UCDv1.2.psl /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain 

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainSort liftovers/ARS-UCDv14_to_ARS-UCDv1.2.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/chainNet liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.2bit.info /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/ars_ucd_14_igc_rmask/ARS-UCD1.0.14.clean.wIGCHaps.fasta.2bit.info liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net /dev/null

/mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/netChainSubset liftovers/ARS-UCDv14_to_ARS-UCDv1.2.net liftovers/ARS-UCDv14_to_ARS-UCDv1.2.sorted.chain liftovers/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain
```

### Liftover of ARS-UCDv14 VCF files to ARS-UCDv1.2

First generate the dictionary of the reference/new assembly or else the picard LiftoverVcf will through a java error that is not very informative as I cannot understand and keep wondering what the problem was!!!

```bash
module load java
module load picard
[kiranmayee.bakshy@assembler2 condensed_vcfs]$ java -jar $PICARD CreateSequenceDictionary R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa O=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa.dict  

# check for the liftover of the smallest chromosome
[kiranmayee.bakshy@assembler2 condensed_vcfs]$ java -jar $PICARD LiftoverVcf I=29.vcf.gz O=liftover_to_v1.2/29.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/ARS-UCDv1.2_to_ARS-UCDv14/ARS-UCDv1.2_to_ARS-UCDv14.liftover.chain REJECT=liftover_to_v1.2/rejected/29.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa WRITE_ORIGINAL_POSITION=true

# Now run a script to liftover all the chr files
[kiranmayee.bakshy@assembler2 condensed_vcfs]$ sbatch -p assemble3 liftover.vcf.sh
```
There seems to be no errors and the script is running...

#### Here is the liftover.vcf.sh script

```bash
#!/usr/bin/sh

module load gatk/3.7  
module load java/jdk1.8.0_92

vcffile=vcf.list
vcflines=`cat $vcffile`

for i in $vcflines
do 
	name=`echo $i | rev | cut -c 8- | rev` 
	echo $name
java -jar /mnt/nfs/nfs1/kiranmayee.bakshy/picard/build/libs/picard.jar LiftoverVcf I=$i O=liftover_to_v1.2/$name.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/ARS-UCDv14_to_ARS-UCDv1.2/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain REJECT=liftover_to_v1.2/rejected/$name.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa WRITE_ORIGINAL_POSITION=true
	echo "java -jar /mnt/nfs/nfs1/kiranmayee.bakshy/picard/build/libs/picard.jar LiftoverVcf I=$i O=liftover_to_v1.2/$name.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/ARS-UCDv14_to_ARS-UCDv1.2/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain REJECT=liftover_to_v1.2/rejected/$name.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa WRITE_ORIGINAL_POSITION=true"
done
```

VCF Liftover was 95% successful. I have concatenated the individual vcf files to make one large vcf.

```bash
module load bcftools
module load samtools
module load htslib

bcftools concat -a -f vcf.list | bcftools annotate -I 'bakshy_%CHROM\_%POS\_%REF\_%FIRST_ALT' |  bgzip  > combined.vcf.gz

tabix -C combined.vcf.gz

bcftools query -l combined.vcf.gz > sample.list
```

Now I need to liftover the haplotype segment coordinates bed file to the current assembly coordinates to start filtering...

```bash
[kiranmayee.bakshy@assembler2 condensed_vcfs]$ /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver seg_cord.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/ARS-UCDv14_to_ARS-UCDv1.2/ARS-UCDv14_to_ARS-UCDv1.2.liftover.chain liftover_to_v1.2/seg.cord.mapped liftover_to_v1.2/seg.unmapped
Reading liftover chains
Mapping coordinates
```

Out of 587 segments, 555 (94.5%) have lifted over to the current version, 32 segments failed to liftover because they were all partially deleted in new.

Now start the filtration...

```bash
#First select only biallelic SNP and get rid of INDels and multiallelic SNPs in the regions of interest (seg.cord.mapped)...
bcftools view -m2 -M2 -v snps combined.vcf.gz -R seg.cord.mapped -O z -o filtration/f1.vcf.gz

#Index for further filtration
tabix -C filtration/f1.vcf.gz

#Remove duplicate entries
/home/kiranmayee.bakshy/vcflib/bin/vcfuniq filtration/f1.vcf.gz | bgzip > filtration/f2.vcf.gz 
tabix -C filtration/f2.vcf.gz

#Now run the pipeline for progressive filtration and generate stats

```

This is the [link](https://github.com/bkiranmayee/My_Labnotes/blob/master/CDDR/Filtration_of_ARS-UCDv1.2_coordinates.md) to the variant filtration document...





