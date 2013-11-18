#!/bin/bash

# Set up Bash anv Vi
# Note: JAVA_HOME is defined in bash_profile. Might need to update path.
cp $HOME/miseq_pipe/software/profiles/bash_profile $HOME/.bash_profile
cp $HOME/miseq_pipe/software/profiles/virc $HOME/.virc

# Go into software and unpack tools
cd $HOME/miseq_pipe/software

# Clean existing unzips
rm -rf $HOME/miseq_pipe/software/jasperstarter
rm -rf $HOME/miseq_pipe/software/samtools-0.1.19
rm -rf $HOME/miseq_pipe/software/bedtools-2.17.0
rm -rf $HOME/miseq_pipe/software/setuptools-1.1.6

tar -jxf samtools-0.1.19.tar.bz2
tar -jxf jasperstarter-2.0.0-bin.tar.bz2
tar -zxf BEDTools.v2.17.0.tar.gz
tar -zxf setuptools-1.1.6.tar.gz
tar -zxf vcftools_0.1.7.tar.gz

# Set up jasperstarter
cd $HOME/miseq_pipe/software/jasperstarter/bin
#mv jasperstarter jasperstarter.old
#mv jasperstarter.exe jasperstarter
chmod a+x jasperstarter.exe

# Set up samtools
cd $HOME/miseq_pipe/software/samtools-0.1.19
cp $HOME/miseq_pipe/software/samtools.Makefile.fix $HOME/miseq_pipe/software/samtools-0.1.19/Makefile
make clean
make
mkdir -p bin
cp samtools.exe bin/samtools

# Set up vcftools (filterVCF.sh)
cd $HOME/miseq_pipe/software/vcftools_0.1.7
make clean
make
cd $HOME/miseq_pipe/software/vcftools_0.1.7/bin
cp vcftools.exe vcftools

# Set up bedtools (coverageBed)
cd $HOME/miseq_pipe/software/bedtools-2.17.0
make clean
make

# Set up python's setuptptools
cd $HOME/miseq_pipe/software/setuptools-1.1.6
python setup.py install

# Set up Python Modules
cd $HOME/miseq_pipe/software/

python get-pip.py
pip install wheel

pip install --use-wheel --no-index --find-links=$HOME/miseq_pipe/software/python_wheel beautifulsoup4
pip install --use-wheel --no-index --find-links=$HOME//miseq_pipe/software/python_wheel lxml
pip install --use-wheel --no-index --find-links=$HOME/miseq_pipe/software/python_wheel illuminate

# Set up windows wrapper
mkdir -p /usr/local/bin
cd $HOME/miseq_pipe/software/windows_wrapper/
cp bin/* /usr/local/bin
cp *.lnk /cygdrive/c/Users/$USER/Desktop
