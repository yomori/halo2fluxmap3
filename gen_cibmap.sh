#!/bin/bash

for i in {0..0}
#for i in {0,1,2,3,4,5,6,7,17,19,21,22,23,26,27}
do
#python convertcatalog.py ${i}
cp params_template.txt input.txt
echo inputfile  = \'catalogs/halocatalog_buzzard-${i}.pksc\'  >> input.txt
echo outbase    = \'/project2/chihway/sims/buzzard/cib/buzzardall\' >> input.txt
python buzzard2marcelo.py
python makemaps.py input.txt
done





