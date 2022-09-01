#!/bin/sh

cat bed_list.txt | xargs -P 5 -I{} sh -c 'findMotifsGenome.pl ../../output_bed/CA/{}.bed hg38 {} -size 200 -mask -nomotif'
