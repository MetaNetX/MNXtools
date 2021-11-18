# A minimalistic but helpful test suite

To execute all tests declared in this directory simply run

	`ALLTESTS.sh`

Principles:

* Test are stored in the `*.test` files.
* Every test consists of one or several UNIX commands. 
* New test must start with a TEST=<test-name> UNIX statement
* A test is succesful when the exit code is zero for every command.

The following UNIX environment variables are used to write tests:

	TEST		name of the test (in the . Setting a new name
    	        starts a new test. Test name must be unique
				in the overall collection

	REF			reference dir (where test input and output data
			    are stored under subversion)
	TMP			scratch dir

	INIT        If not empty, its content will be executed at
            	the sartup of every test.
	SOFT        root of the metanetx software
	NCORE		number

For example:

cd $HOME/devel/metanetx/tests/compute
export INIT='module add R/3.6.1;'
export NCORE=8
export SOFT=../..
export REF=model/bigg_e_coli_core
export TMP=/scratch/local/weekly/test

./run_test.pl fba.test

To update the reference output of most tests, declare

	TMP=$REF

before running the test.


