#!/bin/sh

set -ex

autoreconf -isv

set +x
echo
echo "------------------------"
echo "Initialized build system"
echo "------------------------"
echo

if [ "x$1" != "x" ]; then
    echo "+ running './configure $*' for you"
    echo
    exec ./configure "$@"
else
    echo "run './configure && make' to compile"
fi
