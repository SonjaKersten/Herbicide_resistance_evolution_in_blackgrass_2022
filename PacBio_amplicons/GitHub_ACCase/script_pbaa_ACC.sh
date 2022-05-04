# USAGE:  # ./script_pbaa_ACC.sh <REF> <acc> <GENE> <accuracy>
# example # ./script_pbaa_ACC.sh $PWD/ref/ACCase_bait.fa DE01321_01 ACCase 0.99

# parse parameteres
REF=$1
acc=$2
GENE=$3
accuracy=$4

# declare other variables
WORKTMP=$PWD
ref=$(echo "${REF}" | rev | cut -d'/' -f1 | cut -c 4- | rev)
Q=q"$(printf "%.0f\n" $(echo "-10*(l(1-$accuracy)/l(10))" | bc -l))"
q=$(printf "%.0f\n" $(echo "-10*(l(1-$accuracy)/l(10))" | bc -l))

inputTMP=$WORKTMP/output_fastq

outputTMP=$WORKTMP/output_pbaa_"$ref"_$Q

# create necessary directories
cd $WORKTMP

mkdir -p $outputTMP

cd $outputTMP 

# activate virtual environment
source activate /ebio/abt6_projects8/alopecurus_genome/bin/anaconda2/envs/pbaa_v1.0.0  

# Other necessary softwares:
# samtools

# start
echo ""
echo "Start script"
date

##################################################
################ pbaa Clustering #################
##################################################

# gunzipping input fastq
rsync -av $inputTMP/$acc*.rep*.$GENE.$Q.fastq.gz $outputTMP
sleep 2
gunzip -f $outputTMP/$acc.rep*.$GENE.$Q.fastq.gz
gunzip -f $outputTMP/"$acc"_r.rep*.$GENE.$Q.fastq.gz
sleep 2

# I opted for concatenating the available FASTQ files 
cat $outputTMP/$acc*.rep*.$GENE.$Q.fastq > $outputTMP/$acc.$GENE.$Q.fastq

echo "Indexing fastq files"
samtools faidx $outputTMP/$acc.$GENE.$Q.fastq

# Create list-of-files
ls $outputTMP/$acc.$GENE.$Q.fastq > $outputTMP/$acc.$GENE.$Q.fofn

echo "Clustering with pbAA..."
pbaa cluster \
	$REF \
	$outputTMP/$acc.$GENE.$Q.fofn \
	--min-read-qv $q \
	$outputTMP/$acc.$Q.$GENE

# Remove large and unneccesary files
rm $outputTMP/$acc*.rep*.$GENE.$Q.fastq*
rm $outputTMP/$acc.$GENE.$Q.fofn
rm $outputTMP/$acc.$GENE.$Q.fastq
rm $outputTMP/$acc.$Q.*_ecr.fasta

conda deactivate

##################################################
################# post-EDITING ###################
##################################################

echo "Collecting metadata..."
grep ">" $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta | sed "s/:/ /g" > $outputTMP/$acc.$Q.$GENE.passed.txt

file1=$outputTMP/$acc.$Q.$GENE.passed.txt

num_of_lines=$(cat $file1 | wc -l)
clusters=($(sed 's/.*cluster-//' $file1 | cut -d'_' -f1))
reads=($(sed 's/.*_ReadCount-//' $file1 | cut -d' ' -f1))
frequency=($(sed 's/.*cluster_freq//' $file1 | cut -d' ' -f2))
diversity=($(sed 's/.*diversity//' $file1 | cut -d' ' -f2))
quality=($(sed 's/.*avg_quality//' $file1 | cut -d' ' -f2))


echo "Selecting right clusters..."
rm $outputTMP/$acc.$Q.$GENE.h*.c*.r*.f*_.fasta
# 1 cluster
if [ $num_of_lines -eq 1 ] ;then
	if [ $reads -ge 25 ] ;then
		grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta > $outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.f${frequency[0]}_.fasta
		grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta > $outputTMP/$acc.$Q.$GENE.h2.c${clusters[0]}.r${reads[0]}.f${frequency[0]}_.fasta
	else
		echo "FAILED.1-clusters.minreads" >  $outputTMP/$acc.$Q.$GENE.FAILED.1-clusters.minreads.err
	fi
# 2 clusters	
elif [ $num_of_lines -eq 2 ]; then
	difference=$(echo "${frequency[0]} ${frequency[1]}" | awk '{print $1-$2}')
	sum_reads=$(IFS=+; echo "$((${reads[*]}))")

	if [ $sum_reads -ge 25 ] && [ $(echo "$difference<=0.50" | bc) -eq 1 ] ;then
		echo "passed 2-clusters difference filter"
		grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta > $outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.f${frequency[0]}_.fasta
		grep -A1 "_cluster-${clusters[1]}" $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta > $outputTMP/$acc.$Q.$GENE.h2.c${clusters[1]}.r${reads[1]}.f${frequency[1]}_.fasta
	else
		echo "FAILED.2-clusters.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.2-clusters.quality.err
	fi
# 0 clusters
elif [ $num_of_lines -eq 0 ]; then

	echo "FAILED.0-clusters.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.0-clusters.quality.err
# 3 or more clusters
else 
	echo "FAILED.3-clusters.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.3-clusters.quality.err
fi


echo "Editing fasta files..."
rm $outputTMP/$acc.$GENE.h*.fasta*

# Edit header
if [ -f "$outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.f${frequency[0]}_.fasta" ]; then
	echo "I found a file"
	sed "s/>.*/>$acc.h1/" $outputTMP/$acc.$Q.$GENE.h1.c*.r*.f*_.fasta > $outputTMP/$acc.$GENE.h1.fasta
	sed "s/>.*/>$acc.h2/" $outputTMP/$acc.$Q.$GENE.h2.c*.r*.f*_.fasta > $outputTMP/$acc.$GENE.h2.fasta

	ls $acc.$Q.$GENE.h*.c*.r*.f*_.fasta > $outputTMP/$acc.$Q.$GENE.all.txt 
fi

# Complementary reverse when needed
rc=`grep -c "^GAA" $outputTMP/$acc.$GENE.h1.fasta`
if [ $rc -eq 1 ]; then
       samtools faidx --reverse-complement --mark-strand no $outputTMP/$acc.$GENE.h1.fasta $acc.h1 > $outputTMP/$acc.$GENE.h1.fastaTMP
       mv $outputTMP/$acc.$GENE.h1.fastaTMP $outputTMP/$acc.$GENE.h1.fasta
fi

rc=`grep -c "^GAA" $outputTMP/$acc.$GENE.h2.fasta`
if [ $rc -eq 1 ]; then
       samtools faidx --reverse-complement --mark-strand no $outputTMP/$acc.$GENE.h2.fasta $acc.h2 > $outputTMP/$acc.$GENE.h2.fastaTMP
       mv $outputTMP/$acc.$GENE.h2.fastaTMP $outputTMP/$acc.$GENE.h2.fasta
fi

# Remove unneccessary files
rm $file1

# End
date
echo "End script"
echo ""



