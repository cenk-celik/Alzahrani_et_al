#!/bin/bash

for fileName in data/IRE1_*_{1..6};

do

samp=`basename ${fileName}`

echo "Processing sample ${samp}"

salmon quant -i index -l A \
	-1 ${fileName}/${samp}_1.fq.gz \
	-2 ${fileName}/${samp}_2.fq.gz \
	-p 8 --validateMappings \
	--gcBias -o quants/${samp}_quant
done
