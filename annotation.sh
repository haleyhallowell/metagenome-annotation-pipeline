#!/bin/sh

##moving to the correct place
cd /scratch/aubhah
##the files used in the study had almost identical file names, with the exception of an individual sample number
r=20

##This pipeline has a built in .txt file that reports at each step whether or not the program completed 
##successfully, in this particular instance it is called 'C4LJourney.txt'. This can be named whatever you want
echo "reads processing, here comes sample $r!" >> C4LJourney.txt

##This was performed on a remote super computer, so here I am calling the progrmam
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastx/0.0.14

### Check the quality of our file
fastx_quality_stats -i CB14TANXX_s4_1_GSLv3-7_${r}_*.fastq  -o before_quality_stats_sample${r}_1.txt
fastx_quality_stats -i CB14TANXX_s4_2_GSLv3-7_${r}_*.fastq  -o before_quality_stats_sample${r}_2.txt

##############################  TRIMMOMATIC #########################
# -b from anywhere in the sequence, if overlap 5' end trim it, if in middle or overlapping 3' end then trims and removes everything after.
# This is the index adapter as written 5' to 3';
# -e is the error rate. 2% so if >50 bp match can have 1 mismatch;
# -m is the length so if reads less than 30bp after trimming then removed
#MAKE SURE YOU HAVE THE TRUSEQ3-PE.fa loaded as a file in the directory. it can ba found here: https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa

##loading up program
module load anaconda
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda_3-4.2.0
#Denote where the program is
TRI_PATH=/opt/asn/apps/anaconda_3-4.2.0_cent/share/trimmomatic/trimmomatic.jar

#moving to the correct directory
cd /home/aubhah/pigproject/analysis/2_Trimmomatic

##CB14TANXX_1 and _2 are the forward and reverse reads, respectively. These will be replaced with your forward and reverse file names
#  run Trimmomatic
java -jar $TRI_PATH PE -phred33 /home/aubhah/pigproject/data/CB14TANXX_s4_1_GSLv3-7_${r}_*.fastq /home/aubhah/pigproject/data/CB14TANXX_s4_2_GSLv3-7_${r}_*.fastq /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_1_pe.fq /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_1_se.fq /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_2_pe.fq /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_2_se.fq ILLUMINACLIP:/home/aubhah/pigproject/analysis/2_Trimmomatic/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20

##checklin line counts to make sure trimmed files are equal
wc -l  /home/aubhah/pigproject/data/CB14TANXX_s4_1_GSLv3-7_${r}_*.fastq >> /home/aubhah/pigproject/analysis/2_Trimmomatic/${r}T_clipped_record.txt;
wc -l  /home/aubhah/pigproject/data/CB14TANXX_s4_2_GSLv3-7_${r}_*.fastq >> /home/aubhah/pigproject/analysis/2_Trimmomatic/${r}T_clipped_record.txt;
wc -l /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_1_pe.fq >> /home/aubhah/pigproject/analysis/2_Trimmomatic/${r}T_clipped_record.txt;
wc -l /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_2_pe.fq >> /home/aubhah/pigproject/analysis/2_Trimmomatic/${r}T_clipped_record.txt;

##Checking Trimomatric output

echo "Checking Trimmomatic run:" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt

a= `/home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_1_pe.fq  | wc -l`
b= `/home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_2_.pe.fq | wc -l`

if [ "$a" -ne "$b" ];
	then
		echo "ERROR: Cleaned paired files are assymetrical. Re-run Trimmomatric" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
	else
		echo "Paired files are symetrical, you're doing amazing sweetie!" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
fi

##Check that our adapters were removed, and quality is better than original samples

fastqc  /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_1_pe.fq -o /home/aubhah/pigproject/analysis/3_FastQC
fastqc  /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_2_pe.fq -o /home/aubhah/pigproject/analysis/3_FastQC

##copying Trimmomatic outputs to correct directory
scp /home/aubhah/pigproject/analysis/2_Trimmomatic/sample${r}_*_pe.fq /home/aubhah/pigproject/analysis/4_RemoveHostSequences

##moving to correct directory for next step
cd /home/aubhah/pigproject/analysis/4_RemoveHostSequences

######## Get rid of host sequences #######################
#Alright, so now that we've gotten through the easy part, now we are going to filter out the host sequences to clean our stuff up more
#We will do this by mapping to the pig genome and then keeping the sequences that DID NOT map! 

##loading programs
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load samtools/1.2

#For pig genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCA_000003025.6_Sscrofa11.1/GCA_000003025.6_Sscrofa11.1_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/*.fna.gz
gunzip *.fna.gz
cat *.fna > Sus_scrofarefgenome.fa
bwa index -p Sus_scrofarefgenome.fa -a bwtsw Sus_scrofarefgenome.fa
samtools faidx Sus_scrofarefgenome.fa

##Mapping (using 'bwa mem' for its increased accuracy)
#MAPPING PAIRED END READS
bwa mem -t 4 -P Sus_scrofarefgenome.fa sample${r}_1_pe.fq sample${r}_2_pe.fq > /home/aubhah/pigproject/analysis/4_RemoveHostSequences/sample${r}_pe.remove.sam ;


###Check things worked

echo "Checking BWA run:" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt
if [ -s sample${r}_pe.remove.sam ]; #TRUE=File exists and is bigger than zero
	then
		echo "Merged paired .sam files successfully created." >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
	else
		echo "ERROR: failed to create .sam files, rerun BWA" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
		exit
fi

###Okay so now that we've mapped to the host lets pull out the sequences that DIDNT map
#load the SAMtools environment

#moving to correct directory
cd /home/aubhah/pigproject/analysis/5_PrepforAssembly

samtools view -@ 4 -buh -f12 sample${r}_pe.remove.sam > /home/aubhah/pigproject/analysis/5_PrepforAssembly/sample${r}_pe.unmapped.bam ;
rm sample${r}_pe.remove.sam

##Convert to fastq
 samtools bam2fq sample${r}_pe.unmapped.bam > sample${r}_pe.unmapped.fq 
rm sample${r}_pe.unmapped.bam

#Pull apart
paste - - - - - - - - < sample${r}_pe.unmapped.fq \
    | tee >(cut -f 1-4 | tr "\t" "\n" > sample${r}_1_rfa.fq) \
    |       cut -f 5-8 | tr "\t" "\n" > sample${r}_2_rfa.fq
    
##Checking output
 echo "Checking Samtools:" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt
 if [ -s sample${r}_2_rfa.fq ]; #TRUE=File existis and is bigger than zero
 	then
 echo "ERROR: something went wrong with in this samtools mess" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
 	exit
 fi

#moving files to correct directory
scp sample${r}_pe.unmapped.fq /home/aubhah/pigproject/analysis/6B_IDBAAssembly
rm sample${r}_pe.unmapped.fq

cd /home/aubhah/pigproject/analysis/6B_IDBAAssembly


#################Assembly with IDBA-UD######################
#Loading programs
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2018.12
module load fastx/0.0.14

#convert to fasta 
fastq_to_fasta -v -n -i sample${r}_pe.unmapped.fq -o sample${r}_pe.unmapped.fa
rm sample${r}_pe.unmapped.fq
/opt/asn/apps/anaconda_3-2018.12/bin/idba fq2fa /home/aubhah/pigproject/analysis/5_PrepforAssembly/sample${r}_pe.unmapped.fq sample${r}_pe.unmapped.fa
/opt/asn/apps/anaconda_3-2018.12/bin/idba idba_ud -r sample${r}_pe.unmapped.fa -o sample${r}


########Map back to assembly with Bowtie2########
#Let's get to mapping! First step is to build your reference, second step is to map back to your assembly
bowtie2-build /home/aubhah/pigproject/analysis/6B_IDBAAssembly/sample${r}/contig.fa IDBAsample${r}_contig
bowtie2 -x IDBAsample${r}_contig -1 sample${r}_1_rfa.fq -2 sample${r}_2_rfa.fq -S IDBAbowtie${r}.sam
rm sample${r}_*_rfa.fq 

##output moved to correct directory
scp /home/aubhah/pigproject/analysis/8_Bowtie/IDBAbowtie${r}.sam /home/aubhah/pigproject/analysis/9_Metaphlan
rm sample${r}_*_rfa.fq

#move to correct directory
cd /home/aubhah/pigproject/analysis/9_Metaphlan

##You made it! Time to annotate with Metaphlan3!
samtools view -S -b IDBAbowtie${r}.sam > bowtie${r}.bam
rm IDBAbowtie${r}.sam
samtools bam2fq bowtie${r}.bam > bowtie${r}.fq
rm bowtie${r}.bam
/opt/asn/apps/anaconda_3-2020.07/bin/metaphlan bowtie${r}.fq --input_type fastq --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /home/aubhah/mpa_v20_m20/mpa_v30_CHOCOPhlAn_201901 -t rel_ab_w_read_stats  --ignore_archaea --ignore_eukaryotes --nproc 4 --o bowtie21.fq.bowtie2out.txt_bacteriaout.txt

echo "Checking Metaphlan2 run:" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt
if [ -a /home/aubhah/pigproject/analysis/9_Metaphlan/sample${r}metaphlanout.txt]; #TRUE=File existis 
	then
		echo "Bish your sample is finished. you're amazing!!! I love you!!" >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
	else
		echo "ERROR: metaphlan2 didn't work.. go back and check your code." >> /home/aubhah/pigproject/logs/C4LJourney${r}.txt;
		exit
fi


























