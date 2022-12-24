# PAgonist

<br>
---
title: "<center> RClass </center> "
author: 
  - name: "<center> Sung Young Kim, MD, PhD^1^</center>"
    affiliation: "<center> ^1^Department of Biochemistry, School of Medicine, Konkuk University</center> "
output: 
    html_document:
        toc: false  
        toc_float: 
           collapsed: false
        toc_depth: 3
        mathjax: https://d3eoax9i5htok0.cloudfront.net/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML
---
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

> - <font size='2'> highlight content 1 <b>[$y=ax+b+c$](#here)</b></font>  
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




<div align="center">
```{r, echo = FALSE, fig.width=8, fig.cap='Fig. 1', warning = FALSE , message = FALSE, collapse=TRUE }
source("~/MS/F2.R")
awViz <- function(proj.title = "test", fontname='Arial', fontsize=14, nodewidth= 4, penwidth=1, arrowsize=1, node.color="darkslateblue", text.color="white", customInput = list(letters[1:9], letters[11:19], c(1:9)), awInput = list(letters[1:9], letters[11:19], c(1:9)), seed=1234, aw=ranseq(awInput, "random", seed=seed), awCol=as.list(sample( .col("workflow")(6), 6))) {
require(DiagrammeR);
gvvar  <- list()
gvvar['basic.par'][[1]] = c(fontname=fontname, fontsize=fontsize, nodewidth= nodewidth, penwidth=penwidth, arrowsize=arrowsize, node.color=node.color, text.color=text.color)
gvvar['customInput'][[1]] = customInput
gvvar['awInput'][[1]] = awInput
gvvar['awCol'][[1]] = awCol
gvvar <<- gvvar
  
graph <- grViz( diagram <- "digraph example {
graph [fontname=@@1-1, nodesep = .3, ranksep = .5, rankdir=LR, splines=ortho, compound=true, label='', labelloc=t, labeljust=left, fontsize=30, fontcolor=darkslateblue, layout = dot] # concentrate=false, fontname='Helvetica Bold', layout = dot
node [fontname = @@1-1, fontcolor=@@1-7, shape = rectangle, style=filled, color=@@1-6, fontsize = @@1-2, width=@@1-3, height = .7, margin = '0.3, 0.16', fixedsize=TRUE] 
edge [fontname = @@1-1,  fontcolor=@@1-7, color = '#00000090', dir = both, arrowtail = none, arrowsize = @@1-5, fontsize=@@1-2, penwidth=@@1-4, style='dashed']
                  
############## EDIT ====                
subgraph cluster0 {color='#00000010'; style='filled'; fillcolor= '#00000010'; penwidth=@@1-4; margin=0; label='Pathway Dysregulation Score'; labelloc= 't'; labeljust=l; fontname=@@1-1; fontsize=14; fontcolor='#000000';
sub0[margin = '0, 0', width=1, fixedsize=FALSE, label=<<table width='100%' bgcolor='#00000030' cellborder='0' cellspacing='0' cellpadding='10'> 
<TR> <TD BALIGN='LEFT' BGCOLOR='#000000'> Pathifier <br/>PathTacerIPS___________ </TD><TD BALIGN='LEFT' BGCOLOR='#00000030'> Pathifier<br/>PathTacerIPS  </TD></TR> 
<TR> <TD BALIGN='LEFT' BGCOLOR='#00000030'> Pathifier <br/>PathTacerIPS___________ </TD><TD BALIGN='LEFT' BGCOLOR='#000000'> Pathifier<br/>PathTacerIPS  </TD></TR>  
<TR> <TD BALIGN='LEFT' BGCOLOR='#000000'> Pathifier <br/>PathTacerIPS___________ </TD><TD BALIGN='LEFT' BGCOLOR='#00000030'> Pathifier<br/>PathTacerIPS  </TD></TR>  
</table>> ]   }

'Meta-analysis'[width=1.5,shape=circle,peripheries=4] 'Meta-analysis'->'ML module'
'Pathway-modules'[width=1.5,shape=circle,peripheries=4] 'Pathway-modules'-> sub0:nw->'ML module':nw [style=vis, lhead = cluster0, ltail = cluster0, weight=1]
                  
{rank=same;'Meta-analysis' -> 'Pathway construction (unsupervised)' ->'Pathway-modules'   'Meta-analysis' -> 'Pathway construction (supervised)' ->'Pathway-modules'} 
{rank=same; a[color='#FFFFFF', fontcolor=black, height=.3,label='Features (genes)'] a->'Meta-analysis'} 
{rank=same;'Evaluation'->'Additional Analysis'}
{rank=same;'ML module'-> 'Evaluation'[constraint=true]  }
                  
 ############## EDIT END ==== 
                  
# a[label='@@1-3'] p[shape=point, width=0] a->p[dir=none] p->b[style='dashed'] p->ESR
# gg0 [group=e1, label = 'gg0 label']  # gg0 -> {gg1 gg2}   # ub3l:se-> End:w [weight=2]
# gg1 [label = 'gg1 label\\l <b>.....</b>\\r', fillcolor = pink ]
# edge [color = grey, style=vis, minlen=2, constraint = true ] # invis
# { rank = same; gg1 -> gg2 [label='edge labe bra']} # node align : 1.rank: horizonal > Group (vertical로 노드에 group=e2 등으로...) > weight (edge에 weight=2...) 
}
[1]: gvvar['basic.par'][[1]] # list도 상관없음.
[2]: gvvar['customInput'][[1]]
[3]: gvvar['awInput'][[1]]
[4]: gvvar['awCol'][[1]]
")
graph
   }
awViz()

```
</div>


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
