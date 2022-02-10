# Useful reading: http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

export TMP=/tmp/MNXtools_test/$USER/tests/$$ # scratch dir
mkdir -p $TMP
chmod a+wrx $TMP
#NOTE to use binaries and cache from the container:
# export DEFAULT_PATH=/usr/local/MNXtools
# export OS_PATH='docker run --name mnxtools --mount type=bind,source=/tmp/MNXtools_test,target=/tmp/MNXtools_test --rm -i -t sibswiss/mnxtools'
export MNX_OS_PATH=''
export MNX_DEFAULT_PATH=..

for PERL in perl;
do
    for MNET in \
        ${DEFAULT_PATH:=$MNX_DEFAULT_PATH}/examples/ECOLI/bigg_iML1515 ; # ../examples/HUMAN/metatlas_HumanGEM
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

