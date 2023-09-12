#!/usr/bin/env bash

LD_LIBRARY_PATH=${HOME}/.local/lib
CASADI_PATH=$(python -c "import casadi; print(casadi.__path__[0])")

sudo cp ${CASADI_PATH}/libcasadi.so.3.7 ${LD_LIBRARY_PATH}
sudo cp ${CASADI_PATH}/libcasadi.so ${LD_LIBRARY_PATH}
