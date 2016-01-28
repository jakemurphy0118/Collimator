#!/bin/bash

for i in {1..2}
do
bsub ./CALIOPE $i
done
