#!/bin/bash
USAGE="Usage: $0 11"
if [ $# -ne 1 ]
then
    echo $USAGE
    exit
fi

if [ ! -d 'File' ]
then
    echo "Creating File"
    mkdir 'File'
fi
sftp "ftp:/star/u/hfeng1/Kochratio/UrQMD/ReadUrQMD/File/$1.root" "File/"
