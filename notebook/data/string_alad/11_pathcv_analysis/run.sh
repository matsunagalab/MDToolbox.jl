#!/bin/bash

for i in $(seq 1 16); do sed -e 's/<num>/'${i}'/' pathcv.in >pathcv.${i}.in ; done

for i in $(seq 1 16); do /Users/yasu/genesis/bin/pathcv_analysis pathcv.${i}.in; done

for i in $(seq 1 16); do   awk '{print $1, $3}' ${i}.pathcv >${i}.pathdist; done

