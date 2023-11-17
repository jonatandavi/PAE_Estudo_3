#!/bin/bash

for d in {1..64};
do
    div=$(echo "$d % 3" | bc)
    if [ $div = "0" ]; then
        wait
    fi    
    Rscript scriptMod.R $d &
done