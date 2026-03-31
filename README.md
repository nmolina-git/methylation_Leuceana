# DupGenFinder Pipeline for *Leucaena trichandra*

This pipeline prepares *Leucaena trichandra* genome annotation and protein files for use with DupGenFinder. It also includes preparation of an outgroup species, *Medicago truncatula*, and running DIAMOND searches required for duplication mode analysis.

## Requirements

Before beginning, make sure the following programs are installed:

- DIAMOND
- Perl
- wget
- gzip
- awk
- sed
- DupGenFinder

## 1. Set Up Working Directory

Create a project directory and move into it:

```bash
mkdir DupGenFinder_Project
cd DupGenFinder_Project
```

Create subdirectories for input files and results:

```bash
mkdir input_files
mkdir results
```

Move your *Leucaena trichandra* files into the `input_files` directory:

```bash
mv Ltri1.0_gene_models_Noelle_gendup.gff3 input_files/
mv Ltri1.0_prot_Noelle_gendup.fasta input_files/
cd input_files
```

## Output

Results from DupGenFinder will be saved in the `results` directory inside your working directory.

## Notes

- Make sure gene names in FASTA and GFF files match exactly.
- Check DIAMOND version compatibility.
- Consider using `nohup` for long-running DIAMOND jobs.

## 2. Prepare *Leucaena trichandra* Genome Files

Rename the genome FASTA and GFF files to match the format required by DupGenFinder:

```bash
mv Ltri1.0_gene_models_Noelle_gendup.gff3 Ltr.gff
mv Ltri1.0_prot_Noelle_gendup.fasta Ltr.fasta
```

Create a DIAMOND reference database:

```bash
diamond makedb --in Ltr.fasta -d Ltr
```

Run DIAMOND blastp against the reference database:

```bash
nohup diamond blastp -d Ltr -q Ltr.fasta -o Ltr.blast -p 12 --sensitive --max-target-seqs 5 --evalue 1e-10 --quiet
```

## 3. Prepare Outgroup (*Medicago truncatula*)

Download the outgroup protein FASTA file:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa.gz
gzip -d GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa.gz
mv GCF_003473485.1_MtrunA17r5.0-ANR_protein.faa Mtr.fasta
```

Download the outgroup GFF file:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff.gz
gzip -d GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff.gz
mv GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff Mtr.gff
```

## 4. Format Outgroup Files

The FASTA and GFF files must have matching gene names and headers.

Format the GFF file:

```bash
zgrep -v '^#' GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff | \
awk -F'\t' '$3=="mRNA"{match($9,/Name=([^;]+)/,a); if(a[1]!="") print $1,$3,$4,$5,a[1]}' OFS="\t" > Mtr.gff

awk '!/^#/ {OFS="\t"; print $1, $5, $3, $4}' Mtr.gff > Mtr.gff_edited

rm Mtr.gff
mv Mtr.gff_edited Mtr.gff
```

Format the FASTA file:

```bash
sed -E 's/^>(\S+).*/>\1/' Mtr.fasta > Mtr.fasta_edited

rm Mtr.fasta
mv Mtr.fasta_edited Mtr.fasta
```

## 5. Combine GFF Files and Run DIAMOND

Combine the *Leucaena* and *Medicago* GFF files:

```bash
cat Ltr.gff Mtr.gff > Ltr_Mtr.gff
```

Create a DIAMOND database for the outgroup:

```bash
diamond makedb --in Mtr.fasta -d Mtr
```

Run DIAMOND blastp between *Leucaena trichandra* and *Medicago truncatula*:

```bash
nohup diamond blastp -d Mtr -q Ltr.fasta -o Ltr_Mtr.blast -p 12 --sensitive --max-target-seqs 5 --evalue 1e-10 --quiet
```

## 6. Move Extra Files

Move unnecessary intermediate files into a separate directory:

```bash
mkdir extra_files
mv Ltr.dmnd Mtr.dmnd extra_files/
mv Ltr.fasta Mtr.fasta Mtr.gff extra_files/
```

## 7. Run DupGenFinder

Run DupGenFinder with the prepared input files:

```bash
perl /Programs/DupGen_finder/DupGen_finder.pl \
  -i /media/Scratch8TB/CLG2023/molinan/Leucaena_trichandra_data/Dup_gene_stuff/input_files/ \
  -t Ltr \
  -c Mtr \
  -o results
```

Results will be saved in the `results` directory.
