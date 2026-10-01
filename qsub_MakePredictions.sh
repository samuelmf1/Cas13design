#!/bin/bash
#$ -cwd

DIR=$1

TargetGenes=$(ls ./${DIR}/*fasta)

for target in ${TargetGenes}; do
	echo "Rscript ./scripts/RfxCas13d_GuideScoring.R  ${target} ./data/Cas13designGuidePredictorInput.csv true" | qsub -V -cwd -e /dev/null -o /dev/null
done

