#!/usr/bin/env bash

cd ../bin/BRB-seqTools/releases/

InputDir=${1}
OutputDir=${1}/demultiplexed
mkdir -p ${OutputDir}

#BarcodeFile=../../../data/fastq/barcodes/barcodes_96_V4C_brb.txt
BarcodeFile=${2}
UMILength=${3}


ls ${InputDir}/*_R1.fastq.gz

for Library in ${InputDir}/*_R1.fastq.gz
do
	LibraryFile=${Library##*/} # "file.txt"
	LibraryName=${LibraryFile%_R1.fastq.gz} # "file"
	
	mkdir -p ${OutputDir}/${LibraryName}

	echo "${InputDir}/${LibraryName}_R1.fastq.gz"
	echo "${InputDir}/${LibraryName}_R2.fastq.gz"

	java -jar BRBseqTools-1.6.1.jar Demultiplex \
		-r1 ${InputDir}/${LibraryName}_R1.fastq.gz -r2 ${InputDir}/${LibraryName}_R2.fastq.gz \
		-c ${BarcodeFile} -p BU -UMI ${UMILength} \
		-o ${OutputDir}/${LibraryName}
done

# To run: 
# ./demultiplex_brbseq_fastq.sh ../../../data/fastq/AG0012 ../../../data/fastq/barcodes/barcodes_v4_set0_brb.txt 13
# ./demultiplex_brbseq_fastq.sh ../../../data/fastq/AMP0020 ../../../data/fastq/barcodes/barcodes_96_V4C_brb.txt 9
# ./demultiplex_brbseq_fastq.sh ../../../data/fastq/AMP0027 ../../../data/fastq/barcodes/barcodes_96_V4C_brb.txt 9