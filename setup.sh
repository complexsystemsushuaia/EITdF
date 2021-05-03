#!/bin/bash
#
# Setup script for EITdF
#
# Daniel Badagnani
# Instituto de Ciencias Polares, Ambiente y Recursos Naturales
# Universidad Nacional de Tierra del Fuego
# Ushuaia, Argentina
#
################################################################

if [[ -d "../runEITdF" ]]
then
  echo
  echo "directory ../runEITdF already present"
  echo
else
  mkdir ../runEITdF
fi
cp demo.draft ../runEITdF/demo.ipynb

echo "*****************************************"
echo "**                                     **"
echo "**    EITdF is ready for running       **"
echo "**                                     **"
echo "**      Please use ../runEITdF         **"
echo "**                                     **"
echo "*****************************************"
