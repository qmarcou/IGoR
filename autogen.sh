#!/bin/bash

#First need to build the configure script for jemalloc using autogen
cd libs/jemalloc
bash autogen.sh
cd ../..
CURRENT_DIR=`pwd`
echo Back in $CURRENT_DIR directory

#Rename jemalloc configure.ac file to prevent it from messing with autoconf configuration
if [ -e "./libs/jemalloc/configure.ac" ]
then
echo "Temporarily Renaming jemalloc configure.ac to prevent it from messing with autoconf"
mv ./libs/jemalloc/configure.ac ./libs/jemalloc/configure.ac_tmp
fi

#Run the autotools in order to obtain configure scripts
autoreconf --install

if [ -e "./libs/jemalloc/configure.ac_tmp" ]
then
mv ./libs/jemalloc/configure.ac_tmp ./libs/jemalloc/configure.ac
fi


#Create the markdown file from asciidoc
#https://tinyapps.org/blog/nix/201701240700_convert_asciidoc_to_markdown.html
asciidoc -b docbook README.adoc
pandoc -f docbook -t markdown_strict README.xml -o README.md

# Create a man page and html README from the markdown manual using Pandoc
#adapted from :https://stackoverflow.com/questions/7599447/less-style-markdown-viewer-for-unix-systems#7603703
pandoc -s -f markdown -t html README.md -o docs/README.html
pandoc -s -f markdown -t man README.md -o igor.1

# Create the Doxygen documentation
doxygen doxygen.config
cd ./docs/latex/
make all
cd ../..
mv ./docs/latex/refman.pdf ./docs/IGoR_CPP_manual.pdf

