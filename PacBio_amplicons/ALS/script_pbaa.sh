# USAGE:  # ./script_pbaa.sh <REF> <acc> <GENE> <accuracy>
# example # ./script_pbaa.sh $PWD/ref/ALS_bait.fa DE01321_01 ALS 0.99

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

# activate virtual environment
source activate pbaa_v1.0.0
# Other necessary softwares:
# samtools

# start
echo ""
echo "start script"
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

# Remove unneccesary files
rm $outputTMP/$acc*.rep*.$GENE.$Q.fastq*
rm $outputTMP/$acc.$GENE.$Q.fofn
rm $outputTMP/$acc.$GENE.$Q.fastq
rm $outputTMP/$acc.$Q.*_ecr.fasta


##################################################
################# post-EDITING ###################
##################################################

echo "Concatenating and cleaning outputs..."
cat $outputTMP/$acc.$Q."$GENE"_passed_cluster_sequences.fasta \
	$outputTMP/$acc.$Q."$GENE"_failed_cluster_sequences.fasta \
	> $outputTMP/$acc.$Q.$GENE.TMP1.fasta

grep -B1 "^NNN" $outputTMP/$acc.$Q.$GENE.TMP1.fasta \
	| diff - $outputTMP/$acc.$Q.$GENE.TMP1.fasta \
	| sed '/^> /!d;s/^> //' \
	> $outputTMP/$acc.$Q.$GENE.TMP2.fasta

grep -A1 "ALS_001361Frc" $outputTMP/$acc.$Q.$GENE.TMP2.fasta \
	| diff - $outputTMP/$acc.$Q.$GENE.TMP2.fasta \
	| sed '/^> /!d;s/^> //' \
	> $outputTMP/$acc.$Q.$GENE.all.fasta

rm $outputTMP/$acc.$Q.$GENE.TMP*.fasta

# Produce and edit output test
grep ">" $outputTMP/$acc.$Q.$GENE.all.fasta | sed "s/:/ /g" > $outputTMP/$acc.$Q.$GENE.all.txt

sed -i "s/486280-489800_cluster-/cluster /g" $outputTMP/$acc.$Q.$GENE.all.txt
sed -i "s/_ReadCount-/ ReadCount /g" $outputTMP/$acc.$Q.$GENE.all.txt

echo "Select appropriate clusters..."

# Collect metadata
file1=$outputTMP/$acc.$Q.$GENE.all.txt

num_of_lines=$(cat $file1 | wc -l)
clusters=($(awk '{print $3}' $file1))
reads=($(awk '{print $5}' $file1))
frequency=($(awk '{print $13}' $file1))
diversity=($(awk '{print $15}' $file1))
quality=($(awk '{print $17}' $file1))

# 1 cluster
if [ $num_of_lines -eq 1 ] ;then
	if [ $reads -ge 25 ] && [ $(echo "$frequency>=0.98" | bc) -eq 1 ] && [ $(echo "$diversity<=0.4" | bc) -eq 1 ] && [ $(echo "$quality>=0.7" | bc) -eq 1 ] ;then

		grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.fasta
		grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h2.c${clusters[0]}.r${reads[0]}.fasta
	else
		echo "FAILED.c1.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.c1.quality.err
	fi

# 2 clusters 
elif [ $num_of_lines -eq 2 ]; then
	difference=$(echo "${frequency[0]} ${frequency[1]}" | awk '{print $1-$2}')
	sum_reads=$(IFS=+; echo "$((${reads[*]}))")

	if [ $sum_reads -ge 25 ] && [ $(echo "$difference<=0.85" | bc) -eq 1 ] ;then

		if [ $(echo "${diversity[0]}<=0.4" | bc) -eq 1 ] && [ $(echo "${diversity[1]}<=0.4" | bc) -eq 1 ] && [ $(echo "${quality[0]}>=0.7" | bc) -eq 1 ] && [ $(echo "${quality[1]}>=0.7" | bc) -eq 1 ] ;then

			grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.fasta
			grep -A1 "_cluster-${clusters[1]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h2.c${clusters[1]}.r${reads[1]}.fasta
		else
			echo "FAILED.c2.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.c2.quality.err
		fi
	else
		echo "FAILED.c2.difference" > $outputTMP/$acc.$Q.$GENE.FAILED.c2.difference.err
	fi

# 3 clusters
elif [ $num_of_lines -ge 3 ]; then
	addition=$(echo "${frequency[0]} ${frequency[1]}" | awk '{print $1+$2}')
	difference=$(echo "${frequency[0]} ${frequency[1]}" | awk '{print $1-$2}')
	sum_reads=$(IFS=+; echo "$((${reads[*]}))")

	if [ $(echo "$addition>= 0.96" | bc) -eq 1 ] && [ $(echo "$difference<=0.85" | bc) -eq 1 ];then

		if [ $sum_reads -ge 25 ] && [ $(echo "${diversity[0]}<=0.4" | bc) -eq 1 ] && [ $(echo "${diversity[1]}<=0.4" | bc) -eq 1 ] && [ $(echo "${quality[0]}>=0.7" | bc) -eq 1 ] && [ $(echo "${quality[1]}>=0.7" | bc) -eq 1 ] ;then

			grep -A1 "_cluster-${clusters[0]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.fasta 
			grep -A1 "_cluster-${clusters[1]}" $outputTMP/$acc.$Q.$GENE.all.fasta > $outputTMP/$acc.$Q.$GENE.h2.c${clusters[1]}.r${reads[1]}.fasta
		else
			echo "FAILED.c3.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.c3.quality.err
		fi
	else
		echo "FAILED.c3.addition" > $outputTMP/$acc.$Q.$GENE.FAILED.c3.addition.err
	fi

# 0 clusters
elif [ $num_of_lines -eq 0 ]; then
	echo "FAILED.c0.quality" > $outputTMP/$acc.$Q.$GENE.FAILED.c0.quality.err
# 3 or more clusters
else
	echo "FAILED.cN.somethingelse" > $outputTMP/$acc.$Q.$GENE.FAILED.cN.somethingelse.err
fi

echo "Editing fasta files..."

# Edit header
if [ -f "$outputTMP/$acc.$Q.$GENE.h1.c${clusters[0]}.r${reads[0]}.fasta" ]; then
	echo "I found a file"
	sed "s/>.*/>$acc.h1/" $outputTMP/$acc.$Q.$GENE.h1.c*.r*.fasta > $outputTMP/$acc.$GENE.h1.fasta
	sed "s/>.*/>$acc.h2/" $outputTMP/$acc.$Q.$GENE.h2.c*.r*.fasta > $outputTMP/$acc.$GENE.h2.fasta
fi

# Complementary reverse
# Pattern used to be: "GTG", but changed to "GAA" with newer version of pbaa

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



conda deactivate


# End
date
echo "End script"
echo ""



