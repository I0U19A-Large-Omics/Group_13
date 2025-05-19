#!/bin/bash
echo "Running benchmarks..."


# Run make clean to remove any existing output from previous runs
# Run a make rule that generates a file test_time/test_serial_scored_snps.tsv
# Wrap the make command in time to measure execution time
# Redirects both stdout and stderr using 2>&1 | tee serial.time:
# - 2>&1: merges error output with normal output
# - serial.time: This file will contain the output of the time command, including the real time taken to execute the make command.
# - tee: This command reads from standard input and writes to both standard output and the specified file.
echo "1. Serial processing"
make clean
mkdir -p test_time_report
(time make test_time/test_serial_scored_snps.tsv) 2>&1 | tee test_time_report/serial.time

echo "2. 4 parallel processes"
make clean
mkdir -p test_time_report
(time make test_time/test_p4_scored_snps.tsv) 2>&1 | tee test_time_report/p4.time

echo "3. 8 parallel processes"
make clean
mkdir -p test_time_report
(time make test_time/test_p8_scored_snps.tsv) 2>&1 | tee test_time_report/p8.time

echo "4. 16 parallel processes"
make clean 
mkdir -p test_time_report
(time make test_time/test_p16_scored_snps.tsv) 2>&1 | tee test_time_report/p16.time

echo "5. 36 parallel processes"
make clean
mkdir -p test_time_report
(time make test_time/test_p36_scored_snps.tsv) 2>&1 | tee test_time_report/p36.time
