#!/bin/sh

module load python
virtualenv ve
. ve/bin/activate

for line in $(cat requirements.txt)
do
  pip install $line
done
