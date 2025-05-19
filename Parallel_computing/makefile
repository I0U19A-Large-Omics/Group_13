# In a Makefile, rules are stated in the form of:
# target (output): prerequisites(input)
#	commands 
# The commands are executed when the target is out of date with respect to the prerequisites.
# The target is out of date if it does not exist or if any of the prerequisites are newer than the target.
# The commands are executed in the order they are listed, 
# and the target is considered up to date when all commands have been executed successfully.

# In a Makefile, $@ refers to the target of the rule, while $< refers to the first prerequisite (input file) of the rule.

#Import and copy the input annotated VCF file to the current working directory
data/snps.annotated.vcf: /lustre1/project/stg_00079/teaching/data/snps.annotated.vcf
	mkdir -p data
	cp $< $@ 

# Extract only the SNPs from the annotated VCF file (getting rid of lines with ##)
data/snps_only.tsv: data/snps.annotated.vcf
	mkdir -p data
	bash extract_snps.sh $< $@ 

# Only keep the first 36 lines of the SNPs file for testing purposes
test/test_snps_only.tsv: data/snps_only.tsv
	mkdir -p test
	head -n 36 $< > $@

# Create a header file for the VCF file
data/headers.txt: data/snps.annotated.vcf
	mkdir -p data
	grep "^#" $< > $@

# Uses fake_enformer_score_snps.sh to score the SNPs in the snps_only.tsv file
# cat $< : outputs the file line by line, i.e., each SNP in the file.
# xargs -L 1 : takes each line of the input file and passes it as an argument to the command that follows.
# -P 12 : runs 100 processes in parallel.
# -I{} : replaces {} with the argument passed to the command.
# bash fake_enformer_score_snps.sh "{}" : runs the script fake_enformer_score_snps.sh with the SNP as an argument.
# The output of the script is redirected to the file $@, which is the target of the rule.
data/scored_snps.tsv: data/snps_only.tsv
	mkdir -p data
	cat $< | xargs -L 1 -P 36 -I{} bash fake_enformer_score_snps.sh "{}" > $@

# Create a complete VCF file with the headers and the scored SNPs, run with make complete_scored_snps.vcf
data/complete_scored_snps.vcf: data/headers.txt data/scored_snps.tsv
	mkdir -p data
	cat $^ > $@

# Uses fake_enformer_score_snps.sh to score the SNPs in the test file
test/test_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test
	cat $< | xargs -L 1 -P 36 -I{} bash fake_enformer_score_snps.sh "{}" > $@

# Create a complete test VCF file with the headers and the scored SNPs, run with make test/test_complete_scored_snps.vcf
test/test_complete_scored_snps.vcf: data/headers.txt test/test_scored_snps.tsv
	mkdir -p test
	cat $^ > $@




# REPORT TIME
# Serial processing (no parallelization), you can run this with make test_time/test_serial_scored_snps.tsv
test_time/test_serial_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test_time
	time (cat $< | while read line; do bash fake_enformer_score_snps.sh "$$line"; done > $@)

# Benchmark different parallelization settings, you can run each with make test_time/test_pX_scored_snps.tsv, where X is the number of parallel processes.
# Of course does not make sense to run with lot of parallel processes if you have only a few SNPs. using "nproc" you see we have 36 cores
test_time/test_p4_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test_time
	time (cat $< | xargs -L 1 -P 4 -I{} bash fake_enformer_score_snps.sh "{}" > $@)

test_time/test_p8_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test_time
	time (cat $< | xargs -L 1 -P 8 -I{} bash fake_enformer_score_snps.sh "{}" > $@)

test_time/test_p16_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test_time
	time (cat $< | xargs -L 1 -P 16 -I{} bash fake_enformer_score_snps.sh "{}" > $@)

test_time/test_p36_scored_snps.tsv: test/test_snps_only.tsv
	mkdir -p test_time
	time (cat $< | xargs -L 1 -P 36 -I{} bash fake_enformer_score_snps.sh "{}" > $@)


# Needed in benchmark
# to run benchmarking, you need to run the following command:
# chmod +x benchmark.sh
#./benchmark.sh

.PHONY: clean
clean:
	rm -f *_scored_snps.tsv *.time