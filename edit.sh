#!/bin/sh

sed 'N; 
s/\(.*\)\n\(.*\)/\2\
\1/' dirs.dat > newdirs.dat