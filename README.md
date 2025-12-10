Here is detailed instructions from original source:
https://github.com/WangLabPKU/MAPIT-seq/blob/main/README.md

Install

git clone https://github.com/WangLabPKU/MAPIT-seq
cd MAPIT-seq
conda env create -f env_specific.yml # or use env.yml, more flexible but slower 
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedGraphToBigWig
chmod +x Mapit src/MAPIT-seq.sh bedGraphToBigWig

Optional: Add Mapit to your conda environment PATH

conda activate Mapit-seq
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d

cat <<EOF > $CONDA_PREFIX/etc/conda/activate.d/activate-mapit.sh
#!/usr/bin/env sh
export PATH="\$PATH:/home/gangx/apps/Mapit-seq"
EOF

cat <<EOF > $CONDA_PREFIX/etc/conda/deactivate.d/deactivate-mapit.sh
#!/usr/bin/env sh
export PATH=\$(echo "\$PATH" | tr ':' '\\n' | grep -v '^/home/gangx/apps/Mapit-seq\$' | paste -sd:)
EOF

chmod +x $CONDA_PREFIX/etc/conda/activate.d/activate-mapit.sh $CONDA_PREFIX/etc/conda/deactivate.d/deactivate-mapit.sh

Install MAPIT dependencies

REDItools2

cd ..
git clone https://github.com/BioinfoUNIBA/REDItools2 
conda install mpi4py -c bioconda -c conda-forge
cd REDItools2
pip install -r requirements.txt

‼️ Modifications required for REDItools2 scripts ‼️

Please apply the following changes before running the pipeline:

Add from functools import reduce at line 17 in src/cineca/parallel_reditools.py
Replace line 573 in src/cineca/parallel_reditools.py with:
keys = list(chromosomes.keys())
In src/cineca/reditools.py, replace sys.maxint (line 817) with sys.maxsize
In src/cineca/reditools.py, replace "w" (line 912) with "wt"

SAILOR and FLARE

cd ..
git clone https://github.com/YeoLab/FLARE
conda install snakemake -c bioconda -c conda-forge
pip install deeptools gffutils pyfaidx Bio

Download Reference Sequence and Annotation

cd "your_ref_path"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.p13.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.annotation.gff3.gz

# RepeatMasker
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz  # for mouse, just replace hg38 to mm10/mm39
gzip -d *


# ERCC spike-in
wget https://assets.thermofisher.cn/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip

dbSNP

genomeVersion=GRCh38
mkdir "your_ref_path"/${genomeVersion}_SNP
cd "your_ref_path"/${genomeVersion}_SNP

# human (hg38/GRCh38)
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz # download data in VCF/GATK fold; column 1 starts with "chr"

# mouse (mm10/GRCm38)
wget https://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz

Exome Variant Server or European Variation Archive

# human (hg38/GRCh38)
wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
tar -xvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
for i in {1..22} X Y; do
awk '{if(substr($1, 1, 1) == "#"){print $0}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0; else print "chr"$0 }}}' ESP6500SI-V2-SSA137.GRCh38-liftover.chr${i}.snps_indels.vcf | gzip > EVS_split_chr/chr${i}.gz
done

‼️ The original database link is currently inaccessible. You may download the dataset using the following command and then unzip it directly:

# human (hg38/GRCh38)
wget -O EVS_split_chr.zip "https://zenodo.org/records/17089899/files/EVS_split_chr.zip?download=1"
unzip EVS_split_chr.zip

# mouse (mm10/GRCm38)
wget http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_3/by_species/mus_musculus/GRCm38.p4/GCA_000001635.6_current_ids.vcf.gz
mkdir EVA_split_chr
zcat GCA_000001635.6_current_ids.vcf.gz | awk -v dir_SNP=EVA_split_chr '{if(substr($1, 1, 1) == "#"){print $0 > "EVA_header"}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0 > dir_SNP"/"$1; else print "chr"$0 > dir_SNP"/chr"$1 }}}'
gzip EVA_split_chr/chr*

# mouse (mm39/GRCm39)
# Not recommended due to the lack of supporting dbSNP data for mm39 in NCBI. Consider lifting over from mm10/GRCm38 if needed.
wget http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_3/by_species/mus_musculus/GRCm39/GCA_000001635.9_current_ids.vcf.gz

1000 Genomes Project

# human  
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
mkdir 1000genomes_split_chr
zcat ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz | awk -v dir_SNP=1000genomes_split_chr '{if(substr($1, 1, 1) == "#"){print $0 > "1000genomes_header"}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0 > dir_SNP"/"$1; else print "chr"$0 > dir_SNP"/chr"$1 }}}'
gzip 1000genomes_split_chr/chr*

# mouse
chroms=($(grep '>' $genome_fasta |sed 's/>//' | awk '{print $1}' | grep 'chr' ))
for chr in ${chroms[@]}; do touch void_split_chr/${chr}; done
gzip void_split_chr/chr*

‼️ Important: Ensure all chromosome names in the *_split_chr directories begin with the "chr" prefix. If not, manually prepend "chr" to maintain naming consistency across datasets.

MAPIT configuration

Mapit config --genomeVersion GRCh38 \
             --genomeFasta "full_path"/GRCh38.p13.genome.fa \
             --ERCC "full_path"/"ERCC.fa" \
             --species human \
             --outpath "full_path" \
             --genomeAnno "full_path"/gencode.v40.chr_patch_hapl_scaff.annotation.gff3 \
             --rmsk "full_path"/rmsk.txt \
             --dbSNP "full_path"/GRCh38_SNP/All_20180418.vcf.gz \
             --1000Genomes "full_path"/GRCh38_SNP/1000genomes_split_chr \
             --EVSEVA "full_path"/GRCh38_SNP/EVS_split_chr \
             --Reditools "full_path"/Directory_of_RediTools2.0_software \
             --FLARE "full_path"/Directory_of_FLARE_software

FLARE configuration

"full_path_to_FLARE"/workflow_FLARE/scripts/generate_regions.py <full/path/to/your/genome/gtf/file> <genome_name>_regions

Homer configuration

perl configureHomer.pl -install hg38

Here is the link to SRA reading and downloads: 
[
](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1167221&o=acc_s%3Aa)
