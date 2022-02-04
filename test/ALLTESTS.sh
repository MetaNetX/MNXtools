# Useful reading: http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

export TMP=/tmp/MNXtools_test/$USER/tests/$$ # scratch dir
mkdir -p $TMP
export MNX_DEFAULT_CACHE=../cache/ChemSpace.bindump

for PERL in perl;
do
    for MNET in \
        ../examples/ECOLI/bigg_iML1515 ; # ../examples/HUMAN/metatlas_HumanGEM
    do
        export MNET;
        export PERL;
        perl ./run_test.pl convert.test
        if [ $? != 0 ]
        then
            echo "command failed";
            exit 1;
        fi
    done
done

