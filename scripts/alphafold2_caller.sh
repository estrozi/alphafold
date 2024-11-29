#!/bin/sh
umask 113;
if [ -z "${IBSJOBNAME}" ]; then
  export IBSJOBNAME=fromtty
fi
sudo /storage/Alphafold/scripts/.alphafold2_callee.csh `id -u` "$IBSJOBNAME" "$@"
