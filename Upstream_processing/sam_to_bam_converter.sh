#sam to bam conversion
#MF 10/31/2018

cd /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles

for sample in *.sam
do
	echo $sample
	describer=$(echo ${sample} | sed 's/.sam//')
	echo $describer
	
	#convert file from SAM to BAM format
	samtools view -bS $sample -o ${describer}.uns.bam

	#Sort BAM file
	samtools sort ${describer}.uns.bam ${describer}

	#Index BAM file
	samtools index ${describer}.bam

	#Remove intermediate files
	rm ${describer}.uns.bam

done

#getting a list of bam files
ls *.bam > Culex_bamfiles.txt

#moving sam files to new directory to avoid confusion
mkdir sam_files

mv *.sam ./sam_files

