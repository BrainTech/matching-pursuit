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
    echo
    echo "INFORMATION: --enable-static does not work on MacOS."
    echo "If you want to statically link FFTW library on MacOS, it is necessary"
    echo "to specify the path to FFTW static library using --with-fftw3."
    echo "No dynamic version of the FFTW library may exist in the same directory."
fi
