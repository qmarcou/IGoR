#!/bin/sh
if [ $# == 1 ]; then
	#Automate release number change in configure.ac
	DOTTEDVERSION=$(echo $1 | sed s/-/./g) 
	MYREGEXP='[0-9+\.]+'
	EMAILREGEXP='([[:alnum:]_.-]+@[[:alnum:]_.-]+?\.[[:alpha:].]{2,6})'
	URLREGEXP='(http|https)\://[a-zA-Z0-9\-\.]+\.[a-zA-Z]{2,3}(/[a-zA-Z0-9\-\.]*)*'
	PREVLINE=$(grep -E 'AC_INIT\(\[igor\]\,\ \['$MYREGEXP'\]\,\ \['$EMAILREGEXP'\]\,\ \[igor\]\,\ \['$URLREGEXP'\]\)' configure.ac)
	NEWLINE=$(grep -E 'AC_INIT\(\[igor\]\,\ \['$MYREGEXP'\]\,\ \['$EMAILREGEXP'\]\,\ \[igor\]\,\ \['$URLREGEXP'\]\)' configure.ac | sed -r s/$MYREGEXP/$DOTTEDVERSION/)
	#Now replace it in configure.ac
	#reescape everythin
	ZZ=$(echo $PREVLINE | sed s/'\['/'\\['/g) #escape opening bracket
	ZZ=$(echo $ZZ | sed s/'\]'/'\\]'/g) #escape closing bracket
	ZZ=$(echo $ZZ | sed s/'\.'/'\\.'/g) #escape backslash
	ZZ=$(echo $ZZ | sed 's/)/\\)/g') #escape closing parenthesis
	ZZ=$(echo $ZZ | sed 's/(/\\(/g') #escape opening parenthesis
	ZZ=$(echo $ZZ | sed 's/\:/\\:/g') #escape colum
	ZZ=$(echo $ZZ | sed 's/\ /\\ /g') #escape spaces
	#ZZ=$(echo $ZZ | sed 's/\//\\//g') #escape forward slashes => removed as changing the delimiters to = in the next line is cleaner

	sed -r "s=$ZZ=$NEWLINE=" configure.ac > tmpfile
	rm configure.ac
	mv tmpfile configure.ac
	
	#Automate release number change in README.md
	sed -r s/'Latest released version: '[0-9\.]+/'Latest released version: '$DOTTEDVERSION/ README.adoc > tmpfile
	rm README.adoc
	mv tmpfile README.adoc

	#Commit
	git add configure.ac
	git add README.adoc
	COMMITMESSAGE="IGoR v"$DOTTEDVERSION" release commit."
	git commit -m "$COMMITMESSAGE"

	#Create the packaged archive
	MYPATH=$(pwd)
	NEWDIRPATH=$MYPATH/../igor_$1
	echo Copying repository as $NEWDIRPATH...
	cp -r $MYPATH $NEWDIRPATH
	echo Removing all git related stuff inside
	rm -rf $NEWDIRPATH/.git*
	cd $NEWDIRPATH
	bash autogen.sh
	./configure
	make distclean
	cd ..
	zip -r -D igor_$1.zip ./igor_$1
	echo Cleaning up temporary directory...
	rm -rf $NEWDIRPATH
	echo IGoR v$DOTTEDVERSION release successfully created!
	
	
elif [ $# == 0 ]; then
	echo No release number has been provided... Leaving without performing any action
	echo Please provide a release number as argument with dash seperated version number e.g 5-2-3
else
	echo Too many arguments were passed
	echo Please provide only the release number as argument
fi
