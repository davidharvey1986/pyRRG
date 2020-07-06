#Run an self-enclosed test of pyRRG
virtualenv -p /usr/local/bin/python3.7 bugFixPyRRG
stilts_binary=`which stilts.sh`
silts_binary_path=`dirname $stilts_binary`
sex_binary=`which sex`
sex_binary_path=`dirname $sex_binary`
cd bugFixPyRRG
source activate
export PATH=$(getconf PATH)
export PATH=$PATH:${PWD}/bugFixPyRRG/bin
export PATH=$PATH:$stilts_binary_path
export PATH=$PATH:$sex_binary_path
pip3.7 install numpy
pip3.7 install matplotlib
pip3.7 install pyRRG
