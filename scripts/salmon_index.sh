#!/bin/bash

salmon index -t gencode.v42.transcripts.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
