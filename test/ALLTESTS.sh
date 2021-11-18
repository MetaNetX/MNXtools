# Useful reading: http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

export TMP=/tmp/$USER/tests/$$ # scratch dir
mkdir -p $TMP

for MNET in \
    ../examples/ECOLI/bigg_iML1515;
do
    export MNET # used by convert.test
    perl ./run_test.pl convert.test
    if [ $? != 0 ]
    then
        echo "command failed";
        exit 1;
    fi
done


