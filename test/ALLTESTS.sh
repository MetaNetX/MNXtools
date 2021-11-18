# http://web.archive.org/web/20110314180918/ http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

export TMP=/tmp/$USER/tests/$$
mkdir -p $TMP
export ROOT=..
## INIT will be executed by run_test.pl
export INIT='module add Cheminformatics/LinearProgramming/gurobi/9.0.3;module add R/3.6.1;export GRB_LICENSE_FILE=../../gurobi/$USER/$HOSTNAME/gurobi.lic'
export NCORE=48
export LC_ALL=C

for REF in \
    $ROOT/examples/ECOLI/bigg_e_ \
    $ROOT/examples/ECOLI/ ;
do
    export REF
    #export TMP=$REF
    ./run_test.pl convert.test 
    if [ $? != 0 ]
    then
        echo "command failed";
        exit 1;
    fi
done


