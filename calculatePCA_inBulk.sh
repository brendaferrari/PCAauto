#!/bin/bash

echo starting PCA analysis;

for dir in input/*; 

do cp {PCAauto/calculatePCA.py,options.py,plot.py} "$dir";

cd "$dir";

shopt -s extglob;

if ls -l *.{pdb,xtc};

then echo Starting "${dir##*/}" PCA analysis;

python calculatePCA.py; 

rm -f {calculatePCA.py,options.py,plot.py};

cd ../..;

echo "${dir##*/}" PCA analysis ended;

else echo No .pdb or .xtc files on directory. Please verify; 

cd ../..; fi; done