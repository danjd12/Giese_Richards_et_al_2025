#!/bin/bash

usage(){
cat <<EOF

usage: $0 -w [wild-type sample] -m [mutant sample] -h

This pipeline is intended to go from FASTQ paired-end reads to annotated files
describing variants in the mutant sample that are absent from the wild-type
control.  This is for "sibling subtraction" to identify variants that came
a mutagenesis screen.

Options:
    -w  Wild-type sample, the siblings that lack the phenotype which is being
        assayed for.
    -m  Mutant sample, the siblings that have the desired phenotype.
    -g  Path to genome.fa file
    -t  Path to GTF file with gene annotations
    -h  Show help message (which is this)

EOF
}

while getopts hw:m: options;
do
    case $options in
        h) usage && exit 1 ;;
        w) WT=$OPTARG ;;
        m) MUT=$OPTARG ;;
        g) GENOME=$OPTARG ;;
        t) GTF=$OPTARG ;;
    esac
done

if [[ -z $WT ]] || [[ -z $MUT ]]
then
    usage && exit 1
fi

# The genomes for WS235+ are all the same, so aligning to WS274 is fine
bwa mem $GENOME <(gzip -dc $WT.R1.trm.fastq.gz) <(gzip -dc $WT.R2.trm.fastq.gz) | samtools view -S -b | samtools sort -n -o $WT.name.bam;
samtools fixmate -m $WT.name.bam $WT.fixmate.bam;
rm $WT.name.bam;
samtools sort -o $WT.withdup.bam $WT.fixmate.bam;
rm $WT.fixmate.bam;
samtools markdup -r $WT.withdup.bam $WT.sorted.bam;
rm $WT.withdup.bam;
samtools index $WT.sorted.bam;

# Use bcftools to call variants, identify all and homozygous variants
bcftools mpileup -p -f $GENOME $WT.sorted.bam | bcftools call -vc -O v > $WT.vcf;
bcftools view -g hom -O v $WT.vcf -o $WT.homozygous.vcf;

# Repeat of the above, but on the mutant sample
bwa mem $GENOME <(gzip -dc $MUT.R1.trm.fastq.gz) <(gzip -dc $MUT.R2.trm.fastq.gz) | samtools view -S -b | samtools sort -n -o $MUT.name.bam;
samtools fixmate -m $MUT.name.bam $MUT.fixmate.bam;
rm $MUT.name.bam;
samtools sort -o $MUT.withdup.bam $MUT.fixmate.bam;
rm $MUT.fixmate.bam;
samtools markdup -r $MUT.withdup.bam $MUT.sorted.bam;
rm $MUT.withdup.bam;
samtools index $MUT.sorted.bam;
bcftools mpileup -p -f $GENOME $MUT.sorted.bam | bcftools call -vc -O v > $MUT.vcf;
bcftools view -g hom -O v $MUT.vcf -o $MUT.homozygous.vcf;

# Identify common variants
bedtools intersect -header -a $WT.vcf -b $MUT.vcf > $WT.$MUT.common.vcf;
# Subtract common variants, keep those unique to the mutant sample
bedtools subtract -header -a $MUT.vcf -b $WT.$MUT.common.vcf > unique.$MUT.vcf;
# Call the homozygous mutants
bcftools view -g hom -O v unique.$MUT.vcf -o unique.$MUT.homozygous.vcf;

# Use snpEFF to annotate the variants and determine impact for the homozygous mutants
java -Xmx4g -jar ~/snpEff/snpEff.jar WBcel235.86 unique.$MUT.homozygous.vcf > ann.unique.$MUT.homozygous.vcf;
cat ann.unique.$MUT.homozygous.vcf | ~/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].BIOTYPE" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].HGVS_P" > ssm.$MUT.hom.summary.table;

# Ditto, but includes mutants not homozygous
java -Xmx4g -jar ~/snpEff/snpEff.jar WBcel235.86 unique.$MUT.vcf > ann.unique.$MUT.vcf;
cat ann.unique.$MUT.vcf | ~/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].BIOTYPE" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].HGVS_P" > ssm.$MUT.all.summary.table;

# Use VarScan to look at copynumber variants (duplications/deletions)
samtools mpileup -q 1 -f $GENOME $WT.sorted.bam $MUT.sorted.bam | java -jar ~/varscan/VarScan.v2.4.4.jar copynumber varScan --mpileup 1;
java -jar ~/varscan/VarScan.v2.4.4.jar copyCaller output.copynumber --output-file $MUT.copynumber.called --output-homdel-file $MUT.copynumber.called.homdel;

# Use Bedtools inersect to assign genes to homozygous deletions
bedtools intersect -a <(cat $GTF | awk -F "\t" '$3=="gene" {print $0}') -b $MUT.copynumber.called.homdel > $MUT.homdel.genes.bed