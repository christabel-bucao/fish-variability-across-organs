# This script only works with my base environment
# conda deactivate

cd ../bin/GO-Figure 

InputDir=${1} 

mkdir -p ${1}/go-figure
OutputDir=${1}/go-figure

for GOTerms in ${InputDir}/*.txt
do
	FileName=${GOTerms##*/}
	echo "${FileName%.txt}"
	./gofigure-mac --input ${GOTerms} --input_type standard-plus \
		--output ${OutputDir} --outfile_appendix "${FileName%.txt}" \
		--size members --colours user --sort_by user-descending \
		--colour_label "Number of conditions" --max_label 20 --ontology bpo 
done

# To run:
#./run_go_figure.sh ../../results/07_GO_Enrichment/terms
