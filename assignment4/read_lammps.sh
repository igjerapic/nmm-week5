#!/bin/bash

log_file="log.lammps"
for dir in $(find . -name "rate_*"); do
    echo "Parsing ${log_file} in ${dir}"
    (cd  ${dir} && python ../../scripts/read_lammpsLog.py)
done