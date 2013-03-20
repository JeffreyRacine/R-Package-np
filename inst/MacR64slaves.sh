#!/bin/sh
if [ $# -lt 4 ]; then
	echo "Need 4 arguments" > err.log
	exit 1
fi

Rscript=$1
R_HOME=$4

if [ ! -r $Rscript ]; then
	echo $Rscript "does not exist or is not readable!" > err.log
	exit 1
fi

if  [ "$3" = "needlog" ]; then
	hn=`hostname -s`
	$R_HOME/bin/R64 --no-init-file --slave --no-save < $1 > $hn.$2.$$.log 2>&1
else
	$R_HOME/bin/R64 --no-init-file --slave --no-save < $1 > /dev/null 2>&1
fi
