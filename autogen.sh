#!/bin/sh

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
