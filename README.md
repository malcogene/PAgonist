
<p align="center"><img src="https://raw.githubusercontent.com/malcogene/Tmp/master/img/PAgonist.png" style="width:150px;"></p> 
<!--... -->

# PAgonist

> **Install Dependencies**

* If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
* Open the R command line interface, and install package dependencies (if they have not been installed yet):
* Dependencies: R (>= 3.0.0), shiny (>= 0.8.0), AnnotationDbi, and ggplot2


```{r eval = FALSE}
if (!requireNamespace("BiocManager", quietly = T)) install.packages("BiocManager")
BiocManager::install("limma")
install.packages(c("ggplot2", "ggrepel"))
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
devtools::install_github("malcogene/PAgonist")
```



----


> **Package Download**

* Package source file can be downloaded from [here](https://raw.githubusercontent.com/malcogene/PAgonist/main/R/pagonist.R)


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
<br>

Our analysis first starts with determining the connectivities between the gene expression profiles of ranked background gene set $b g$ and the target gene set $t g$ by calculating the Kolmogorov-Smirnoff (K-S) statistic score from the hit and miss scores, which represent the cumulative sum of the ranks of genes that do or do not, respectively, appear in both of the two gene sets. Initially, two K-S scores are calculated, $S_{u p}$ for the upregulated genes in the target gene set and $S_{d o w n}$ for the downregulated genes and are later combined to produce the composite connectivity score, $S_{\text {total }}$.

$bg={g_{1}, g_{2} ..., g_{n}}$, the ranked background gene set of the $N$ genes and $r_{bg}=r\left(g_{j}\right)$, of their expression profiles associated with a target phenotype.

$tg=$ the target gene set

$N_{b g}=$ number of the background gene set

$N_{t g}=$ number of the target gene set

$R_{b g, t g}=\left(r_{b g, t g}\right)^{s}$ gene expression profile found in the background genes and the target gene set, where An exponent $s$ to control the weight of $r_{b g, t g}$

$p=$ penalty value determined by the user

$$eCDF_{Hit, i}=\sum_{0}^{i}\left|R_{bg, tg}\right| *\left(\Sigma\left|R_{bg, tg}\right|\right)^{-1}$$ cumulative sum of the ranks of 'hit' genes found in between and leading up to gene i divided by the total sum of the ranks of all 'hit' genes

$$eCDF_{\text {Miss, } i}=\frac{1}{N_{bg}-N_{tg}} *\left(N_{bg}-N_{tg, j}\right)$$ cumulative sum of the number of 'miss' genes found in between and leading up to gene $\mathrm{j}$ divided by the total sum of the number of all 'miss' genes

The ‘hit’ and ‘miss’ scores are then subtracted to calculate $S_{i}$, the total cumulative score for the ranked genes including and leading up to gene $i$.


$$
S_{i}=S_{Hit,i}-S_{Miss, i}=\sum_{0}^{i}\left|R_{bg, tg}\right| *\left(\Sigma \mid R_{bg, tg}\right)^{-1} \frac{1}{N_{bg}-N_{tg}} *\left(N_{bg}-N_{tg, i}\right)
$$


In order to prioritize drugs that concordantly reverse the diseased profile with minimal unwanted consequences, we multiply a discordant constant $D$ to the score. $D$ shows how unified the direction of the drug-driven "hits" is. Ideally, a drug will upregulate genes that are decreased in the diseased state and vice versa, in which case $D$ will equal to 1 and have no effect on the score. Depending on the sign of the maximum score, we adjust the scores by varying the value of $D$. In cases when the maximum score is positive, which indicates that the target gene set contains downregulated genes and the drug induces a tendency towards upregulation, downregulation of genes represent an unwanted byproduct. Considering that overexpressing genes with drugs is especially difficult, the downregulation of genes that are already downregulated by the disease is highly undesirable, which we penalize with a penalty $p$. In cases when the maximum score is negative, which indicates that target gene set contains upregulated genes and the drug induces a tendency towards downregulation,

$$
D=\frac{1}{2}\left(1+\left(\sum R_{bg, tg}\mid\left(R_{bg, t g}>0\right)\right)-\left(\sum R_{bg, tg}\mid\left(R_{bg, tg}<0\right)\right)\left(\Sigma\left|R_{bg, tg}\right|\right)^{-1}\right)
$$



$$
S_{i, a d j}=  \left(\sum_{0}^{i}\left|R_{b g, t g}\right|{ }^{*}\left(\Sigma\left|R_{b g, t g}\right|\right)^{-1}-\frac{1}{N_{b g}-N_{t g}} *\left(N_{b g}-N_{t g, i}\right)\right)^{p} \quad  \text { if } S_{i, \max }>0
$$

$$
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \left(\sum_{0}^{i}\left|R_{b g, t g}\right| *\left(\Sigma\left|R_{b g, t g}\right|\right)^{-1}-\frac{1}{N_{b g}-N_{t g}} *\left(N_{bg}-N_{tg, i}\right)\right)*\left(1-D\right) \quad  \text { if } S_{i, \max }<0
$$

$$
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \left(\sum_{0}^{i}\left|R_{b g, t g}\right|{ }^{*}\left(\Sigma\left|R_{b g, t g}\right|\right)^{-1}-\frac{1}{N_{b g}-N_{t g}} *\left(N_{b g}-N_{t g, i}\right)\right)*\left(1-D\right) \quad  \text { if } S_{i, \max }<0
$$



The total score is obtained by calculating and subtracting the individual scores for the upregulated and downregulated target gene sets, where the score with the absolute maximum value is selected and adjusted by the process above.

$$
S_{\text {total }}=\left[\begin{array}{ll}
\frac{1}{2} *\left(S_{\text {up, adj, } \max }-S_{\text {down, adj, } \max }\right) & \text { if } \operatorname{sign}\left(S_{\text {up, adj, } \max }\right) \neq \operatorname{sign}\left(S_{\text {down,adj, max }}\right) \\
0 & \text { if } \operatorname{sign}\left(S_{\text {up, adj, } \max }\right)=\operatorname{sign}\left(S_{\text {down,adj, max }}\right)
\end{array}\right.
$$

The significance of the total score is evaluated through a permutation test.

$$
H_{0}: S=S_{\text {total}} 
$$

$$
P=\sum_{0}^{N}(|S|>|S_{total}|)*N^{-1}
$$


We further analyze the genes involved in drug's unwanted effects through overrepresentation analysis (ORA), where we calculate the p-value through a hypergeometric distribution. By identifying genes that are overrepresented in this group, we can explore ways to address them, either with another round of enrichment analysis to discover a supplementary drug to eliminate the unwanted effects or a pathway enrichment analysis to investigate how these unwanted genes may function.

$$
P(X \geq x)=1-P(X \leq x-1)=1-\sum_{i=0}^{x-1} \frac{\left(\begin{array}{l}
M \\
i
\end{array}\right)\left(\begin{array}{l}
N-M \\
n-i
\end{array}\right)}{\left(\begin{array}{l}
N \\
n
\end{array}\right)}
$$
