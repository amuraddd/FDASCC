# FDASCC Package in Python with rpy2 Wrapper
--------------------------------------------
### Documentation  
- The setup requires that the use must have R downloaded 
- They must also have packages gtools, splines2, pracma, RSpectra installed in the R environment (this can be included in the install_fdascc func)
- run the set_r_home function by specifying the location where R is installed to set the r_home environ variable for rpy2 to use
- have the FDASCC folder saved on their machine
- run the install_fdascc function with the location of the FDASCC folder to download the R code as a package
-----------------  
### Set up  
1. Install R Version 3.6.3  
    - [Download R for Windows](https://cran.r-project.org/)
    - [Download R for Mac](https://cran.r-project.org/)
    - [Download R for Linux](https://cran.r-project.org/)

The following steps only need to be done once to set up the base R environment:  

Run the code snippet below - prior to running the code make sure R is installed in a location on your computer where admin priveleges are not required and also make sure the FDASCC R package is downloaded on your computer in a location which is easily accessible:

- The FDASCC R Package can be downloaded from: https://github.com/amuraddd/FDASCC/blob/main/FDASCC/FDASCC.zip

```python
from FDASCC_py.setup import set_r_home, install_fdascc
install_fdascc(path_to_R="<path for where R is downloaded on your computer>", 
              path_to_fdascc="<path where the FDASCC R Package is saved on your computer>")
```

Runnning the Python code above will download the required dependencies along with the FDASCC R package into your R environment which can then be sourced in Python.

-----------------
### Run the exmaples in example.R  

Follow the code in `fdascc_examples.ipynb` to run the examples in `examples.R`. The python wrapper utilizes the library `rpy2` which is a python interface for `R`. The documentation for this library can be found at this [link](https://rpy2.github.io/doc/latest/html/introduction.html). A sample notebook `rpy2_functionality.ipynb` is also provided which shows the basic functionality of the `rpy2` library.