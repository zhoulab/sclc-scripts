#!/bin/sh

module load python
virtualenv ve
module purge
module load python/2.7.6
module load gcc
. ve/bin/activate
for line in $(cat requirements.txt)
do
  pip install $line
done
