#!/bin/bash

command -v Rscript >/dev/null 2>&1

if [ $? -eq 0 ]; then
    echo "Rscript command is available"
    
    Rscript FUNEL_control.r 2>FUNEL_control.err 1>FUNEL_control.log
    if [ $? -eq 0 ]; then
        Rscript FUNEL_sample.r 2>FUNEL_sample.err 1>FUNEL_sample.log
        if [ $? -eq 0 ]; then
            echo "Success: sample_file dealed. Have a nice day!"
        else
            echo "Error: FUNEL_sample.r failed, check FUNEL_sample.err for more information"
            exit 1
        fi
    else
        echo "Error: Rscript FUNEL_control.r failed, check FUNEL_control.err for more information"
        exit 1
    fi
else
    echo "Error: Rscript command is not available, please install Rscript first"
    exit 1
fi

# EOF
