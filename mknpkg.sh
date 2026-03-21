#!/bin/sh

maver=$(cat npver)
miver=${maver##*-}
maver=${maver%-*}

if [ "$1" ] 
    then 
    miver=$((miver+1)) 
    echo "$maver"-"$miver" > npver
fi

echo updating DESCRIPTION
sed -i.bak -e 's/Version[:].*/'"Version: $maver-$miver/" -e 's/Date[:].*/'"Date: $(date +%Y-%m-%d)/" DESCRIPTION

rm DESCRIPTION.bak
