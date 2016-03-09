#!/bin/bash

# Takes a filename and cleans the file for direct `load` into matlab. Saves the
# resulting file in the same directory containing the input file, with
# appending ".cleaned" to the filename.

# TODO document that genes are sorted and to replace their names with their
# index for import into matlab. Probably want to keep a lookup table of index
# -> name for later.

dir=$(dirname "$1")
filename=$(basename "$1")
extension="${filename##*.}"
without_extension="${filename%.*}"

sed -e 's/geneNum//g' $1 > $dir/$without_extension".cleaned".$extension
