#!/bin/bash

for i in $(seq 1 16); do sed -e 's/<num>/'${i}'/' crd_conv.in >crd_conv.${i}.in ; done

for i in $(seq 1 16); do /Users/yasu/genesis/bin/crd_convert crd_conv.${i}.in; done

