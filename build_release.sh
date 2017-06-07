#!/bin/sh
if [ $# == 1 ]; then
	#Automate release number change in configure.ac
	DOTTEDVERSION=$(echo $1 | sed s/-/./g) 
	MYREGEXP='[0-9\.]+'
	EMAILREGEXP='([[:alnum:]_.-]+@[[:alnum:]_.-]+?\.[[:alpha:].]{2,6})'
	PREVLINE=$(grep -E 'AC_INIT\(\[igor\]\,\ \['$MYREGEXP'\]\,\ \['$EMAILREGEXP'\]\)' configure.ac)
	NEWLINE=$(grep -E 'AC_INIT\(\[igor\]\,\ \['$MYREGEXP'\]\,\ \['$EMAILREGEXP'\]\)' configure.ac | sed -r s/$MYREGEXP/$DOTTEDVERSION/)
	echo $NEWLINE
	#Now replace it in configure.ac
	#reescape everythin
	ZZ=$(echo $PREVLINE | sed s/'\['/'\\['/g)
	ZZ=$(echo $ZZ | sed s/'\]'/'\\]'/g)
	ZZ=$(echo $ZZ | sed s/'\.'/'\\.'/g)
	ZZ=$(echo $ZZ | sed 's/)/\\)/g')
	ZZ=$(echo $ZZ | sed 's/(/\\(/g')

	sed -r "s/$ZZ/$NEWLINE/" configure.ac > tmpfile
	rm configure.ac
	mv tmpfile configure.ac
	
	#Automate release number change in README.md
	sed -r s/'Latest released version: '[0-9\.]+/'Latest released version: '$DOTTEDVERSION/ README.md > tmpfile
	rm README.md
	mv tmpfile README.md

	#Commit
	git add configure.ac
	git add README.md
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
	echo Please provide a release number as argument
else
	echo Too many arguments were passed
	echo Please provide only the release number as argument
fi
