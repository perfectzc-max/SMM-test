#env 
conda create --name py2env python=2.7
conda activate py2env
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install pysam

#run main.py
