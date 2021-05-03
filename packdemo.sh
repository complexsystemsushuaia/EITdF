#!/bin/bash
#
# Getting ../runEITdF/demo.ipynb ready for commit
#
# Daniel Badagnani
# Instituto de Ciencias Polares, Ambiente y Recursos Naturales
# Universidad Nacional de Tierra del Fuego
# Ushuaia, Argentina
#
################################################################

if [[ -f "../runEITdF/demo.ipynb" ]]
then
  cp ../runEITdF/demo.ipynb demo.draft
  git add demo.draft
  echo
  echo "  Demo file is ready for commit"
  echo
else
  echo "  ERROR: No demo file for adding to git"
fi
