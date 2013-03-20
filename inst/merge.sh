#!/bin/sh

RMPI_VER=0
RMPI_PREF=Rmpi_

for file in ${RMPI_PREF}*.tar.gz
  do
  if test "$file" \> "$RMPI_VER"
      then
      RMPI_VER=${file}
  fi
done

RMPI_VER=${RMPI_VER#${RMPI_PREF}}
RMPI_VER=${RMPI_VER%.tar.gz}

rm -rf Rmpi
tar -xvzf ${RMPI_PREF}${RMPI_VER}.tar.gz


# identical procedure for np

NP_VER=0
NP_PREF=np_

for file in ${NP_PREF}*.tar.gz
  do
  if test "$file" \> "$NP_VER" 
      then
      NP_VER=${file}
  fi
done

NP_VER=${NP_VER#${NP_PREF}}
NP_VER=${NP_VER%.tar.gz}

rm -rf np
tar -xvzf ${NP_PREF}${NP_VER}.tar.gz


# check for destination dir, if not present, create
# if present, do not clobber

DEST_DIR=npRmpi

if ! test -d ${DEST_DIR} 
    then
    mkdir ${DEST_DIR}
fi

# extract both np and Rmpi to the npRmpi directory
# files in the base install directory will clobber one another
# others that cause problems are exceptions that can be dealt
# with individually

STRIPPER=`gtar --help | grep -- "-strip" | sed -e 's/ *\(-[^=]*\).*/\1=1/'`

gtar -C ${DEST_DIR} ${STRIPPER} -xvzf ${RMPI_PREF}${RMPI_VER}.tar.gz 2> filez.out
gtar -C ${DEST_DIR} ${STRIPPER} -xvzf ${NP_PREF}${NP_VER}.tar.gz 2>> filez.out

# filez.out may now be processed to look for overwritten files
# verify that we know which files are being clobbered

CLOBBERED=`awk '{ print $2 }' filez.out | sort | uniq -d`
CLOBBERED=`echo $CLOBBERED`

PREDICTED="DESCRIPTION R/zzz.R cleanup configure src/Makevars.in src/Makevars.win"

rm filez.out

#if ! test "${CLOBBERED}" = "${PREDICTED}"
#    then
#
#    SPATTERN='s/ /\
#/g'
#    echo Warning: the merge script does not know how to handle the following files from 
#    echo np and Rmpi which conflict:
#    echo `echo "${CLOBBERED} ${PREDICTED}" | sed -e "${SPATTERN}" | sort | uniq -u`
#    exit 1
#fi


## there are files which need to be updated to reflect the new package
## structure. deal with them now.

## this section grabs all top-level functions from the Rmpi files and adds them
## to the namespace

echo -n temp.out

cd npRmpi/R

## the following code can possibly be fooled specially crafted strings in the .R files

for file in R*.R
  do
  sed -e 's/[#].*//' -e 's/[{]/{++^^/g' -e 's/[}]/}--^^/g' -e 's/\([a-zA-Z.][a-zA-Z0-9._]*\) *<- *function/??\1??/g' $file | 
  awk 'BEGIN { 
  FS="([{}]|\\^\\^)" 
  count=0
}
{
  for ( x = 1; x <= NF; x++ ) { 
    if ( $x ~ /[+][+]/ ) { count +=1 }
    else if ( $x ~ /[-][-]/ ) { count -=1 }
    else if ( $x ~ /[?][?].*[?][?]/ && count == 0) { print $x }
  }
}' | sed -e 's/.*[?][?]\([a-zA-Z.][a-zA-Z0-9._]*\)[?][?].*/\1/' >> ../../temp.out
done

cd ../..
sort < temp.out > rmpifuns.out

rm temp.out

awk 'BEGIN { FS="\n" }
{ print }
$0 ~ /export[(].*/{ 
  while (( getline < "rmpifuns.out" ) > 0) { print ",",$0 }
  print ", .Last.lib"
  close("rmpifuns.out")
}' npRmpi/NAMESPACE | sed -e "1s/np/npRmpi/" > NAMESPACE.new

rm rmpifuns.out

mv NAMESPACE.new npRmpi/NAMESPACE

## remove the now unnecessary INDEX file
rm npRmpi/INDEX

## create a new description
sed -i.old -e "s/Package:.*/Package: npRmpi/" npRmpi/DESCRIPTION

rm npRmpi/DESCRIPTION.old

# modify zzz.R
cd npRmpi/R

# convert Rmpi sources to unix files..
for file in R*.R
  do
  tr -d '\r' < $file > tmp.d2u; mv tmp.d2u $file
done

sed -e 's/\(Nonparametric[^"]*\)\\n/\1 + Rmpi '"$RMPI_VER"'\\n/' -e 's/\(library.dynam.unload\)[(]"np"/\1("npRmpi"/' zzz.R |
awk 'BEGIN { FS="\n" }
{ print }
END { 
  print ".onLoad <- function (lib, pkg) {"
  print "  library.dynam(\"npRmpi\", pkg, lib)"
  print "  if(!.Call(\"mpi_initialize\",PACKAGE=\"npRmpi\"))"
  print "    stop(\"Cannot start MPI_Init(). Exit\")"
  print "  if (exists(\".Random.seed\") && "
  print "      round(.Random.seed[1]-5,-1) == .Random.seed[1]-5) {"
  print "    rm(.Random.seed, envir=.GlobalEnv)"
  print "  }"
  print "}"
}' > zzz.R.new

mv zzz.R.new zzz.R

echo '*** Rcoll.old   Mon Jul 31 17:36:03 2006
--- Rcoll.R     Tue Aug  1 09:56:01 2006
***************
*** 55,66 ****
      .Call("bin_nchar", x[1],PACKAGE = "Rmpi")
  }
  
! mpi.bcast.cmd <- function (cmd=NULL, rank=0, comm=1){
      if(mpi.comm.rank(comm) == rank){
          cmd <- deparse(substitute(cmd), width.cutoff=500)
      cmd <- paste(cmd, collapse="\"\"/")
      mpi.bcast(x=nchar(cmd), type=1, rank=rank, comm=comm)
!     invisible(mpi.bcast(x=cmd, type=3, rank=rank, comm=comm))
      } 
      else {
          charlen <- mpi.bcast(x=integer(1), type=1, rank=rank, comm=comm)
--- 55,71 ----
      .Call("bin_nchar", x[1],PACKAGE = "Rmpi")
  }
  
! mpi.bcast.cmd <- function (cmd=NULL, rank=0, comm=1, caller.execute = FALSE){
      if(mpi.comm.rank(comm) == rank){
+       if(caller.execute) tcmd <- substitute(cmd)
          cmd <- deparse(substitute(cmd), width.cutoff=500)
      cmd <- paste(cmd, collapse="\"\"/")
      mpi.bcast(x=nchar(cmd), type=1, rank=rank, comm=comm)
!     bcast.out <- (mpi.bcast(x=cmd, type=3, rank=rank, comm=comm))
!     if (caller.execute) 
!       eval(tcmd, envir = parent.frame())
!     else
!       invisible(bcast.out)
      } 
      else {
          charlen <- mpi.bcast(x=integer(1), type=1, rank=rank, comm=comm)' | patch Rcoll.R


if ! test $?
    then
    exit 1
fi

# make sure the correct package name is always provided

LOF=`grep -il "PACKAGE *= *\"\(np\|Rmpi\)\"" *.R`

for file in $LOF
  do
  sed -i.old -e 's/\(PACKAGE *=\) *"[^"]*"/\1 "npRmpi"/' -e 's/\(package *=\) *"[^"]*"/\1 "npRmpi"/' $file
done

rm *.old
cd ..

# the final step is to make a sane build script skeleton
# rm configure configure.ac configure.in

cd src

RHOME=`R RHOME`
RHOME=`echo $RHOME | sed -e 's/[/]/\\\\\//g'`

SFILES=`ls *.c`
SFILES=`echo $SFILES`

sed -e "s/\(R_INCLUDE_DIR=\).*/\1${RHOME}\/include/" -e "s/\(SOURCES=\).*/\1${SFILES}/" ../inst/Makefile.generic > Makefile
