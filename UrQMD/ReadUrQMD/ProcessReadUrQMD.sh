#!/bin/bash
PWG="$HOME/pwg"
USER=`whoami`
if [ ! -L "$PWG" ]; then
    if [ -L "$HOME/data01" ]; then
        PWG="$HOME/data01"
	elif [ -L "/star/data01/pwg/${USER}" ]; then
		PWG="/star/data01/pwg/${USER}"
    else
        echo "No PWG link found. Maybe You Can Set It Mannually."
        exit
    fi
fi
Energy=0
LOST=0
FileNumber='all'
USAGE="Usage: $0 -e Energy [-f FileNumber] [-L Read Lost.list]"
while getopts "b:e:p:f:L" opt; do
    case $opt in
    e) Energy=$OPTARG ;;
    f) FileNumber=$OPTARG ;;
    L) LOST=1 ;;
    \?) echo "Invalid option: -$OPTARG ${USAGE}" >&2; exit 1 ;;
    esac
done
if [ $Energy -eq 0 ]; then
    echo "Energy are required. ${USAGE}" >&2
    exit
fi
if [ $LOST -eq 1 ]; then
    LIST='Lost.list'
else
    LIST="${Energy}.list"
fi

read -p "Energy: ${Energy}; FileNumber: ${FileNumber}. Continue (y/n)?" CONT
case "$CONT" in
y | Y) echo "Go!" ;;
n | N)
    echo "Nope!"
    exit
    ;;
*) echo "Invalid!" ;;
esac

Prefix=$(pwd | sed "s|$HOME/||")
outDir=$PWG/$Prefix/out
logDir=$PWG/$Prefix/log
errDir=$PWG/$Prefix/err
schedDir=$PWG/$Prefix/generate
for dir in $logDir $errDir $schedDir; do
    if [[ ! -d $dir ]]; then
        echo "Creating the directory $dir"
        mkdir -p $dir && ln -s $dir .
    fi
done

star-submit-template -template Analyze.xml -entities outDir=$outDir,logDir=$logDir,errDir=$errDir,schedDir=$schedDir,workDir=$PWD,Energy=$Energy,FileNumber=$FileNumber,LIST=$LIST

