#!/bin/bash
echo "Running benchmarks..."

echo "1. Serial processing"
make clean
time make test_time/test_serial_scored_snps.tsv 2>&1 | tee serial.time

echo "2. 10 parallel processes"
make clean
time make test_time/test_p10_scored_snps.tsv 2>&1 | tee p10.time

echo "3. 25 parallel processes"
make clean
time make test_time/test_p25_scored_snps.tsv 2>&1 | tee p25.time

echo "4. 50 parallel processes"
make clean 
time make test_time/test_p50_scored_snps.tsv 2>&1 | tee p50.time

echo "5. 100 parallel processes"
make clean
time make test_time/test_p100_scored_snps.tsv 2>&1 | tee p100.time
