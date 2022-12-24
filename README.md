# PAgonist


<!--   주석처리  toc:table of content, html은 각 타이틀 >, md는 ##로 하는 것이 보기좋음    
author: "Sung Young Kim, MD, PhD^[mailto:palelamp@gmail.com]</center>"
date: "`r format(Sys.time(), '%d %B %Y')`"

output: 
    md_document:
        variant: markdown_github     # GitHub Flavored Markdown   github md는 mathjax 안됨
-->

<!--
This page provides the supplementary R code and data to reproduce the experiments in the following paper : **Accurate prediction of acquired EGFR TKIs resistance using a pathway-based individualized machine learning approach**  
-->

> - <font size='2'> highlight content 1 <b>[
> $y=ax+b+c$
> ](#here)</b></font>  
\\[ x = {-b \pm \sqrt{b^2-4ac} \over 2a} \\]

$$ h(t) = h_0(t) \times exp(b_1x_1 + b_2x_2 + ... + b_px_p) $$

$y=ax+b+c$


----

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
* Dependencies: R (>= 3.0.0), shiny (>= 0.8.0), WGCNA, igraph, shinyBS, RColorBrewer, Hmisc, psych, RJSONIO, whisker, yaml, pheatmap, preprocessCore, GO.db, AnnotationDbi, impute, and ggplot2


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




<br><br><br>
<div class="container3"> <div class="child">
```{r echo=FALSE, fig.align='center', message=FALSE, warning=FALSE}
require(ReporteRs)
data(mtcars)
.table(mtcars)

```
</div></div><br><br>



<style>
/* CSS */
div.shadow {
  width: 80%;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19);
  text-align: center;
}
.container1 {
  position: relative;
  width:100%;
  overflow: hidden;
  padding-top: 56.25%; /* 16:9 Aspect Ratio */
}
.container2 {
  position: relative;
  width:100%;
  overflow: hidden;
  padding-top: 100%; 
}
.responsive-iframe {
  position: absolute;
  top: 0;
  left: 0;
  bottom: 0;
  right: 0;
  width: 100%;
  height: 100%;
  border: none;
}
.container3 {
  display: flex;
  justify-content: center;
}
.child {
  background-color: #F6F6F6;
  border: 1px solid black;
  margin: 0 auto;
}
</style>



<br><br>
<center><div class="container1"><iframe class="responsive-iframe" src="https://www.youtube-nocookie.com/embed/vPEa0gNlxNI?rel=0&amp;controls=0&amp;showinfo=0" frameborder="0"></iframe></div></center> 
<br><br>
<center><div class="container2" ><iframe class="responsive-iframe" src="https://centromics.org/centro"  frameborder="0"></iframe></div></center>
<br><br>
<br><br>
<center><div class="shadow"><img src="https://unsplash.it/600.jpg?image=251"></div></center> 
<br><br>
