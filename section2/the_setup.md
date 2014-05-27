# The Setup

Download and install R http://cran.r-project.org/ and RStudio http://www.rstudio.com/ if you do not have them already. Rstudio is optional but it is a great tool if you are just starting to learn R.
You will need specific data sets to run the codes in this document. Download the data.zip[URL to come] and extract it to your directory of choice. The folder name should be “data” and your R working directory should be level above the data folder. That means in your R console, when you type “dir(“data”)” you should be able to see the contents of the data folder. You can change your working directory by *setwd()* command and get your current working directory with *getwd()* command in R2. In RStudio, you can click on the top menu and change the location of your working directory via user interface.


## Installing packages
R packages are add-ons to base R that help you achieve additional tasks that are not directly supported by base R. It is by the action of these extra functionality that R excels as a tool for computational genomics. Bioconductor project (http://bioconductor.org/) is a dedicated package repository for computational biology related packages. However main package repository of R, called CRAN, has also computational biology related packages. In addition, R-Forge(http://r-forge.r-project.org/), GitHub(https://github. com/), and googlecode(http://code.google.com) are other locations where R packages might be hosted.
You can install CRAN packages using install.packages(). (# is the comment character in R)

## Gettting help on functions and packages


### if more help needed...


