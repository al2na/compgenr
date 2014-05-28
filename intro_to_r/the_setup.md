# The Setup

Download and install R http://cran.r-project.org/ and RStudio http://www.rstudio.com/ if you do not have them already. Rstudio is optional but it is a great tool if you are just starting to learn R.
You will need specific data sets to run the codes in this document. Download the data.zip[URL to come] and extract it to your directory of choice. The folder name should be “data” and your R working directory should be level above the data folder. That means in your R console, when you type “dir(“data”)” you should be able to see the contents of the data folder. You can change your working directory by *setwd()* command and get your current working directory with *getwd()* command in R2. In RStudio, you can click on the top menu and change the location of your working directory via user interface.


## Installing packages
R packages are add-ons to base R that help you achieve additional tasks that are not directly supported by base R. It is by the action of these extra functionality that R excels as a tool for computational genomics. Bioconductor project (http://bioconductor.org/) is a dedicated package repository for computational biology related packages. However main package repository of R, called CRAN, has also computational biology related packages. In addition, R-Forge(http://r-forge.r-project.org/), GitHub(https://github. com/), and googlecode(http://code.google.com) are other locations where R packages might be hosted.
You can install CRAN packages using install.packages(). (# is the comment character in R)

```r
# install package named 'randomForests' from CRAN
install.packages("randomForests")
```

You can install bioconductor packages with a specific installer script

```r
# get the installer package
source("http://bioconductor.org/biocLite.R")
# install bioconductor package 'rtracklayer'
biocLite("rtracklayer")
```

You can install packages from github using *install_github()* function from **devtools**

```r
library(devtools)
install_github("roxygen")
```

Another way to install packages are from the source.

```r
# download the source file
download.file("http://goo.gl/3pvHYI", destfile = "methylKit_0.5.7.tar.gz")
# install the package from the source file
install.packages("methylKit_0.5.7.tar.gz", repos = NULL, type = "source")
# delete the source file
unlink("methylKit_0.5.7.tar.gz")
```

You  can also update CRAN and Bioconductor packages.

```r
# updating CRAN packages
update.packages()

# updating bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
```


## Installing packages in custom locations
If you will be using R on servers or computing clusters rather than your personal computer it is unlikey that you will have administrator access to install packages. In that case, you can install packges in custom locations by telling R where to look for additional packages. This is done by setting up an .Renviron file in your home directory and add the following line:
```
R_LIBS=~/Rlibs
```

This tells R that “Rlibs” directory at your home directory will be the first choice of locations to look for packages and install packages (The directory name and location is up to you above is just an example). You should go and create that directory now. After that, start a fresh R session and start installing packages. From now on, packages will be installed to your local directory where you have read-write access.

## Getting help on functions and packages
You can get help on functions by help() and help.search() functions. You can list the functions in a package with ls() function


```r
libray(MASS)
ls("package:MASS") # functions in the package
ls() # objects in your R enviroment
# get help on hist() function
?hist
help("hist")
# search the word "hist" in help pages
help.search("hist")
??hist
@
```

### More help needed?
In addition, check package vignettes for help and practical understanding of the functions. All Bionconductor packages have vignettes that walk you thorugh example analysis. Google search will always be helpful as well, there are many blogs and web pages that have posts about R. R-help, Stackoverflow and R-bloggers are usually source of good and reliable information.




