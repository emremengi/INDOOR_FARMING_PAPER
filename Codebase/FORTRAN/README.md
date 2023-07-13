To run the fortran code for model by iteself (not the file with the GA), run the following at the command line:
$ source /opt/intel/oneapi/setvars.sh intel64
$ ifort {filename} -lm -03 -o {object filename}
$ ./{object filename}
