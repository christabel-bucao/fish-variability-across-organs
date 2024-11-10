#!/usr/bin/env bash

MetadataFile=${1}
FastqDir=${2}

# Check if input arguments are provided
if [[ -z ${MetadataFile} || -z ${FastqDir} ]]; then
	echo "Usage ${0} DRE_metadata.tsv /path/to/fastq/dir"
	exit 1
fi

head ${MetadataFile}

# Step 1: Read metadata file and store [ProjectSeqBatch]-[Barcode]-[SampleName] mappings

## Read metadata file line-by-line
### Extract only relevant fields by column name
ExtractFields=$(awk -F'\t' 'NR==1 {
    # Get relevant column indices
    for (i=1; i<=NF; i++) {
        if ($i=="ProjectSeqBatch") ProjectSeqBatchCol=i
        else if ($i=="Barcode") BarcodeCol=i
        else if ($i=="SampleName") SampleNameCol=i
    }
    # Debug: Print column numbers found
    #print "Debug: ProjectSeqBatch column = " ProjectSeqBatchCol "\t Barcode column = " BarcodeCol, "\t SampleName column = " SampleNameCol | "cat 1>&2"

    # Ensure that all required columns were found
    if (!(ProjectSeqBatchCol && BarcodeCol && SampleNameCol)) {
        print "Error: Required columns not found in header."
        exit 1
    }
} NR>1 { 
    # Print only necessary columns
    print $(ProjectSeqBatchCol), $(BarcodeCol), $(SampleNameCol)
}' ${MetadataFile})

# Debug: Print ExtractFields content
echo "Extracted Fields:"
echo "$ExtractFields"

## Arrays to store [ProjectSeqBatch]_[Barcode] as keys and [SampleName] as values
Keys=()
Values=()

while read -r ProjectSeqBatch Barcode SampleName; do
    # Skip if any of the required fields are empty
    [[ -z ${ProjectSeqBatch} || -z ${Barcode} || -z ${SampleName} ]] && continue

    # Create a unique key by combining ProjectSeqBatch and Barcode
    Key="${ProjectSeqBatch}_${Barcode}"

    # Append the key and sample name to the arrays
    Keys+=("$Key")
    Values+=("$SampleName")

    echo "Adding key ${Key} with sample name ${SampleName}"
    echo "Current number of keys: ${#Keys[@]}, Current number of values: ${#Values[@]}"
done <<< "$ExtractFields"

## Check the number of keys and values added
echo "Number of keys: ${#Keys[@]}"
echo "Number of values: ${#Values[@]}"

## Check that arrays were populated
#for i in "${!Keys[@]}"; do
#    echo "Key: '${Keys[i]}', Value: '${Values[i]}'"  # Print each key and corresponding value
#done

## Function to look up a sample name by ProjectSeqBatch and Barcode
lookup_sample_name() {
    local SearchKey=${1}
    local LocalKeys=${2}
    local LocalValues=${3}

    local FoundValue=""

    #echo "Searching for key: ${SearchKey}"

    # Access array lengths with indirect referencing
    local NumKeys=$(eval "echo \${#${LocalKeys}[@]}")
    local NumValues=$(eval "echo \${#${LocalValues}[@]}")

    #echo "Number of entries in LocalKeys: ${NumKeys}"
    #echo "Number of entries in LocalValues: ${NumValues}"    

    # Loop over indices using indirect referencing
    for ((i=0; i<NumKeys; i++)); do
        local CurrentKey=$(eval "echo \${${LocalKeys}[i]}")
        local CurrentValue=$(eval "echo \${${LocalValues}[i]}")

        #echo "Checking Key: '${CurrentKey}', Value: '${CurrentValue}'"  # Debug each key-value pair
        
        if [[ "${CurrentKey}" == "${SearchKey}" ]]; then
            FoundValue="${CurrentValue}"
            #echo "Found: ${FoundValue}"  # Debug: Print the match and exit the loop
            break
        fi
    done

    if [[ -n ${FoundValue} ]]; then
        echo "${FoundValue}"
    else 
        echo "Sample name not found for key: ${SearchKey}" >&2
    fi
}

## Test lookup function
SearchKey="AG0012_S1_A05"
SampleName=$(lookup_sample_name "${SearchKey}" Keys Values)

echo "Sample name for ${SearchKey} is ${SampleName}"

# Step 2: Loop over each fastq.gz file in current directory
## Extract [ProjectSeqBatch] from directory name

Species=$(basename ${FastqDir} | cut -d'_' -f1)
ProjectSeqBatchDir=$(basename ${FastqDir} | cut -d'_' -f2-)

echo "Species: ${Species}"
echo "ProjectSeqBatchDir: ${ProjectSeqBatchDir}"

for FastqFile in ${FastqDir}/*.fastq.gz; do
    # Extract Barcode (filename without extension)
    FastqBarcode=$(basename ${FastqFile} .fastq.gz)
    echo "FastqBarcode: ${FastqBarcode}"

    # Create the unique key for lookup
    Key="${ProjectSeqBatchDir}_${FastqBarcode}"

    echo "Key: ${Key}"

    # Use the lookup_sample_name() function to get the Sample Name
    SampleName=$(lookup_sample_name "${Key}" Keys Values)
    echo "SampleName: ${SampleName}"

    # Check if SampleName exists for the key
    if [[ -n ${SampleName} ]]; then
        # Construct new file name
        NewFileName="${Species}_${SampleName}_${ProjectSeqBatchDir}_${FastqBarcode}.fastq.gz"

        # Rename the file
        mv ${FastqFile} ${FastqDir}/${NewFileName}
        echo "Renamed ${FastqFile} to ${NewFileName}"
    else
        echo "Warning: No matching entry for Barcode ${FastqBarcode} in ProjectSeqBatch ${ProjectSeqBatchDir}"
    fi
done

# To run:

## Test run:
#./rename_fastq_files.sh ../data/sample_metadata/ELU_sample_annotation.tsv ../data/fastq/test/ELU_AG0012_S1

## Actual run:
#./rename_fastq_files.sh ../data/sample_metadata/ELU_sample_annotation.tsv ../data/fastq/AG0012/demultiplexed/ELU_AG0012_S1

#./rename_fastq_files.sh ../data/sample_metadata/DRE_sample_annotation.tsv ../data/fastq/AMP0020/demultiplexed/DRE_AMP0020_S1

#./rename_fastq_files.sh ../data/sample_metadata/LOC_sample_annotation.tsv ../data/fastq/AG0012/demultiplexed/LOC_AG0012_S0
#./rename_fastq_files.sh ../data/sample_metadata/LOC_sample_annotation.tsv ../data/fastq/AMP0027/demultiplexed/LOC_AMP0027_S2
