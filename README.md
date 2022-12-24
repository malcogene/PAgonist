
# PAgonist

<p align="center"><img src="https://raw.githubusercontent.com/malcogene/Tmp/master/img/faviconCentro10.png" style="width:100px;height:100px;"></p> 

> **Main Dataset**

* The main method function R file can be downloaded from [here](http://centromics.org/info/142sup/mainFunctions.R)
* Preprocessed dataset can be downloaded from [here](http://centromics.org/info/142sup/EGFRTKIs_8set.RData)
* The pathways used for model construction can be downloaded from [here](http://centromics.org/info/142sup/p.KEGG.PID.BioCarta.RData)

----

> **Package Download**

* Package source file can be downloaded from [here](http://centromics.org/info/142sup/mainFunctions.R)


----

> **Install Dependencies**

* If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
* Open the R command line interface, and install package dependencies (if they have not been installed yet):
* Dependencies: R (>= 3.0.0), shiny (>= 0.8.0), AnnotationDbi, and ggplot2


```{r eval = FALSE}
if (!requireNamespace("BiocManager", quietly = T)) install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
install.packages(c("ggplot2", "ggrepel", "WGCNA", "igraph", "Hmisc"))
```
* Please, install versions 0.8.0 for shiny. <!--We are working to make the package compatible with the new versions of the packages as soon as possible.--> To install the recommended versions for shiny, just type the following commands on the R command-line:
```{r eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_version("shiny", "0.8.0")
```


----


> **Install package via Github (Recommended)**

To install the latest version of package via Github, run following commands in R:
```{r eval = FALSE}
devtools::install_github("malcogene/RClass")
```



----


> **Install package from the source**

- **Linux/Mac OS**
    - Download the package *.tar.gz.
    - Open a command prompt. then:
```{r eval = FALSE}
install.packages(path_of_the_downloaded_file, repos = NULL, type="source")
```


- **Windows**
    - Start by reviewing the section on Windows packages in [the R Installation and Administration manual](https://cran.r-project.org/doc/manuals/R-admin.html), then carefully follow the instructions from [The Windows toolset appendix](https://cran.r-project.org/doc/manuals/R-admin.html#The-Windows-toolset).
    - Download the package *.tar.gz.
    - Make sure you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed.
    - Open a command prompt. then:
```{r eval = FALSE}
install.packages(path_of_the_downloaded_file, repos = NULL, type="source")
```


----

