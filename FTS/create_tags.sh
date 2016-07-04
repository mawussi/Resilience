#!/bin/bash

find . -name "*.F90"   -print > taglist.txt
find . -name "*.[Ff]"   -print >> taglist.txt
cat taglist.txt | etags -
rm taglist.txt