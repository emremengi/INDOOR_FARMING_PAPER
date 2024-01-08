**Overview** <br />

**Introduction** <br />

**Results** <br />

**Usage** <br />

**Setup** <br />

**Running** <br />
There are two python versions, one C version, and one Fortran version of this project. For results consistent with those in the publication "A digital-twin and rapid optimization framework for optical design of indoor farming systems," we recommend running the Fortran vrsion. <br />

To run Fortran code, be sure to check your system meets the following requirements.
1. Download and install compiler at https://www.scivision.dev/intel-oneapi-fortran-install/ <br />
2. Download install files from here: https://www.intel.com/content/www/us/en/developer/articles/news/free-intel-software-developer-tools.html <br />
3. Install “Base Kit” and “HPC Kit” to get everything needed <br />
4. Specify compiler location by following these instructions: https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/specifying-the-location-of-compiler-components.html <br />
5. In a terminal window, run source /<install-dir>/setvars.sh intel64 e.g., if the install directory is /opt/intel/oneapi/, run: source /opt/intel/oneapi/setvars.sh intel64 <br />
6. (Optional) If you want this command to always run everytime you open the terminal you can add that line to your ~/.bash_profile or ~/.zshrc (whichever file is sourced in your terminal login shell br />

Once you know your system meets the above requirements, follow the steps below to run any Fortran code in the terminal. <br />
Guide to using the compiler from the command line: https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/invoke-the-compiler.html <br /><br />
$ source /opt/intel/oneapi/setvars.sh intel64 <br />
$ ifort {filename} -lm -03 -o {object filename} <br />
$ ./{object filename} <br />

Now that you know how to run any Fortran file, to **run the indoor farming "lightbox" model by itself,** <br />
1. Navigate to "INDOOR_FARMING_PAPER/Codebase/FORTRAN/MODEL/"
2. Identify the file "PURE-RAY-LIGHTBOX.f" with the latest version number e.g. V13 is later than V12. Then, where <XX> is the latest version number, run <br />
3. ifort PURE-RAY-LIGHTBOX-V<XX>.f -lm -03 -o PURE-RAY-LIGHTBOX-V<XX>.o <br />
4. ./PURE-RAY-LIGHTBOX-V<XX>.o <br />
5. After the executable file runs, a file "PATH.dat" will be created. <br />
6. Using your visualization software of choice (we use Tecplot to produce the plots in the associated paper), you can visualize the data in PATH.dat. <br />

In PATH.dat, each line represents the data in a time step. The columns represent the following data: <br />
1. Column 1: **TODO** <br />
2. Column 2: **TODO** <br />
3. Column 3: **TODO** <br />
4. Column 4: **TODO** <br />
5. Column 5: **TODO** <br />
6. Column 6: **TODO** <br />
7. Column 7: **TODO** <br />
8. Column 8: **TODO** <br />
9. Column 9: **TODO** <br />
10. Column 10: **TODO** <br />
11. Column 11: **TODO** <br />
12. Column 11: **TODO** <br />
12. Column 12: **TODO** <br />

If you desire to edit the model file ("PURE-RAY-LIGHTBOX.f"), we only recommend editing in a few locations:

Now that you know how to run the indoor farming "lightbox" model by itself, to **run the optimization of the indoor farm optical parameters,** <br />
1. Navigate to "INDOOR_FARMING_PAPER/Codebase/FORTRAN/GA/"
2. Identify the file "MACHINE-LEARNING-PURE-RAY-LIGHTBOX.f" with the latest version number e.g. V13 is later than V12. Then, where <XX> is the latest version number, run <br />
3. ifort MACHINE-LEARNING-PURE-RAY-LIGHTBOX-V<XX>.f -lm -03 -o MACHINE-LEARNING-PURE-RAY-LIGHTBOX-V<XX>.o <br />
4. ./MACHINE-LEARNING-PURE-RAY-LIGHTBOX-V<XX>.o <br />
5. After the executable file runs, 
  



**Authors** <br />
Emre Mengi <br />
Carla Becker <br />
Mostafa Sedky <br />
Shao-Yi Yu <br />
Tarek Zohdi <br />

**Contact** <br />
emre_mengi@berkeley.edu <br />
carlabecker@berkeley.edu <br />
msedky@berkeley.edu <br />
syyu410@berkeley.edu <br />
zohdi@berkeley.edu <br />

**Citation** <br />
Mengi, E., Becker, C.J., Sedky M., Yu S., and Zohdi, T.I. A digital-twin and rapid optimization framework for optical design of indoor farming systems.  Computational Mechanics (2023). https://doi.org/10.1007/s00466-023-02421-9 <br />

**License** <br />
No licensing at this time -- all open source. <br />

**Acknowledgments** <br />
Firstly, we would like to thank our advisor Tarek Zohdi for his contributions to this work and for his support and guidance while working on this project. We would also like to acknowledge Jim Pantaleo, the Industry Ambassador at AIFS, and Hanna Bartram, the former Education and Public Engagement Coordinator at AIFS. <br />

**Funding** <br />
This work has been partially supported by the UC Berkeley College of Engineering and the USDA AI Institute for Next Generation Food Systems (AIFS), USDA award number 2020-67021-32855.
