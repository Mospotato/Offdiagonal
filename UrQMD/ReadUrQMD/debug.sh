#!/bin/sh
if [ $# -ne 1 ]; then
    echo "Usage: $0 Energy"
    exit 1
fi
time ./run -f test.list -n DEBUG -e $1
