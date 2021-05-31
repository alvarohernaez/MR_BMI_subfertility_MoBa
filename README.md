# DMRScan
R package for detection of differentially methylated regions with adjustment for multiple testing

# Important Info
## Authors
CM Page, L Vos, BK Andreasen

## Package structure
The function `DMRScan()` requires three key inputs; 
  - a vector of test statistics for each CpG site, 
  - the corresponding positions (chromosome and bp position). 

Additional inputs are the 
  - maximum allowed distance within a region
  - the minimum number of CpGs within a region.

### Installation 
For a stable release from Bioconductor use
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("DMRScan")
```
For the developmental version from Github, use 
```R
devtools::install_github("christpa/DMRScan")
```
Requred R-packages for DMRScan are:
* Matrix 
* MASS 
* RcppRoll 
* GenomicRanges
* IRanges
* methods
* mvtnorm
* stats
* parallel

=======
### Usage
Please see the vignette.

# Citation
To cite the DMRScan package in publications please use:

Page, C. M., Vos, L., Rounge, T. B., Harbo, H. F., & Andreassen, B. K. (2018). Assessing genome-wide significance for the detection of differentially methylated regions. *Statistical applications in genetics and molecular biology.* 


A BibTeX entry for LaTeX users is
```BibTeX
@article{DMR_paper,
  title={Assessing genome-wide significance for the detection of differentially methylated regions},
  author={Page, Christian M and Vos, Linda and Rounge, Trine B and Harbo, Hanne F and Andreassen, Bettina K},
  journal={Statistical applications in genetics and molecular biology},
  year={2018},
  publisher={De Gruyter}
}```


```{r}
sessionInfo()
```

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] DMRScan_0.01.02

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.13            mvtnorm_1.0-6           lattice_0.20-35        
 [4] IRanges_2.10.5          RcppRoll_0.2.2          bitops_1.0-6           
 [7] MASS_7.3-47             GenomeInfoDb_1.12.3     grid_3.4.2             
[10] stats4_3.4.2            zlibbioc_1.22.0         XVector_0.16.0         
[13] S4Vectors_0.14.7        Matrix_1.2-11           RCurl_1.95-4.8         
[16] parallel_3.4.2          compiler_3.4.2          BiocGenerics_0.22.1    
[19] GenomicRanges_1.28.6    GenomeInfoDbData_0.99.0
```


End of README
