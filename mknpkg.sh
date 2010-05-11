#!/bin/sh

maver=`cat npver`
miver=${maver##*-}
maver=${maver%-*}

if [ $1 ] 
    then 
    let 'miver+=1' 
    echo $maver-$miver > npver
fi

echo updating DESCRIPTION
sed -i.bak -e 's/Version[:].*/'"Version: $maver-$miver/" -e 's/Date[:].*/'"Date: `date +%Y-%m-%d`/" DESCRIPTION

rm DESCRIPTION.bak

cd R

echo updating zzz.R
sed -i.bak -e 's/version [0-9.][0-9.]*-[0-9][0-9]*/'"version $maver-$miver/" zzz.R 

rm zzz.R.bak

cd ..