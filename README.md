# Gdf6Y Analysis Scripts

## Structure

This repository contains multiple R scripts generating the data presented in *Richter et al.* 
Each analysis is separately stored in it's folder. For mouse and human the exact R environment that was initially used by the authors, can be restored with the `renv` package and the scripts can be re-run. The *N. furzeri* gonads analysis can also be reconstructed with a little more manual effort.

## Re-running the analysis scripts

Require R versions:

|Analysis | R version|
|:--------|---------:|
|`rnaseq-Mmu` | `4.1.3`|
|`rnaseq-Hsa` | `4.1.3`|
|`rnaseq-Nfu-gonads` | `3.1.3`|

### Running the mouse and human analyses

After cloning this repo to your local directory change into the analysis folder you are interested in (for example `rnsaq-Mmu`) and invoke the associated R version according to the table above.

This automatically activates an `renv` script which downloads and installs itself. Allow also the installation of the `BiocManager` package. 

If R is started in that directory the first time, it requires some attention and time as `renv` downloads all the packages to a global library (*cache*) on your system. Once, those are downloaded, the packages are actually only linked to the analysis directory (this is called *installation*).

Please refer to Appendix I to see the expected behavior after starting R the first time inside an analysis folder.

When all packages are installed for that analysis, `renv::status()` should indicate a `y` for `installed` for all packages. If not, see chapter Troubleshooting.

Then, run the R script via sourcing:

```
source(file = "Mmu-DEG-analysis.R")
```

This produces the Supplementary Data file in xlsx format for the particular species.


### Running the *N. furzeri* gonads analysis

This analysis is divided into two scripts, `Nfu-count_reads.R` and `Nfu-DEG-analysis`. The first can **not** be re-run as it requires the original mapping results (i.e., the `bam` files). The latter starts with the raw counts and should be re-runable under original conditions. This requires R `3.1.3` and the following two packages.

```
library(DESeq2)
library(AnnotationDbi)
```

Versions of all packages according to `sessionInfo()` output.

```
R version 3.1.3 (2015-03-09)
Platform: x86_64-unknown-linux-gnu (64-bit)
Running under: SUSE Linux Enterprise Server 15 SP5

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=C       
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices datasets  utils    
[8] methods   base     

other attached packages:
 [1] AnnotationDbi_1.28.2      Biobase_2.26.0           
 [3] DESeq2_1.6.3              RcppArmadillo_0.5.500.2.0
 [5] Rcpp_0.12.0               GenomicRanges_1.18.4     
 [7] GenomeInfoDb_1.2.5        IRanges_2.0.1            
 [9] S4Vectors_0.4.0           BiocGenerics_0.12.1      

loaded via a namespace (and not attached):
 [1] acepack_1.3-3.3     annotate_1.44.0     base64enc_0.1-3    
 [4] BatchJobs_1.6       BBmisc_1.9          BiocParallel_1.0.3 
 [7] brew_1.0-6          checkmate_1.6.2     cluster_2.0.3      
[10] codetools_0.2-14    colorspace_1.2-6    DBI_0.3.1          
[13] digest_0.6.8        fail_1.2            foreach_1.4.2      
[16] foreign_0.8-66      Formula_1.2-1       genefilter_1.48.1  
[19] geneplotter_1.44.0  ggplot2_1.0.1       grid_3.1.3         
[22] gridExtra_0.9.1     gtable_0.1.2        Hmisc_3.16-0       
[25] iterators_1.0.7     lattice_0.20-33     latticeExtra_0.6-26
[28] locfit_1.5-9.1      magrittr_1.5        MASS_7.3-45        
[31] munsell_0.4.2       nnet_7.3-11         plyr_1.8.3         
[34] proto_0.3-10        RColorBrewer_1.1-2  reshape2_1.4.1     
[37] rpart_4.1-10        RSQLite_1.0.0       scales_0.3.0       
[40] sendmailR_1.2-1     splines_3.1.3       stringi_0.5-5      
[43] stringr_1.0.0       survival_2.38-3     tools_3.1.3        
[46] XML_3.98-1.3        xtable_1.7-4        XVector_0.6.0      
```

When the environment is ready, run the R script via sourcing:

```
source(file = "Nfu-DEG-analysis.R")
```

This script produces three csv files, containing the DESeq results for the three time points.

### Troubleshooting for `renv`

*  In case not all the packages are installed (but in cache/recorded) use the `renv::restore()` command.

*  The message `The project is out-of-sync -- use ``renv::status()`` for details.` after starting R in the analysis directory
can be ignored when `renv::status()` confirms that all packages are installed.

*  For further information, please refer to the manual pages of [renv](https://rstudio.github.io/renv/articles/renv.html).

## Appendix

### I - 1st time behavior of renv

```
The following package(s) are missing entries in the cache:
- annotate
- AnnotationDbi
...
...
...

```

```
These packages will need to be reinstalled.

# Downloading packages -------------------------------------------------------
- Downloading BiocManager from CRAN ...         OK [file is up to date]
Successfully downloaded 1 package in 0.31 seconds.

The following package(s) will be installed:
- BiocManager [1.30.22]
These packages will be installed into "~/R/gdf6y-analysis-scripts/rnaseq-Mmu/renv/library/R-4.1/x86_64-pc-linux-gnu".

Do you want to proceed? [Y/n]: y
```

```
# Installing packages --------------------------------------------------------
- Installing BiocManager ...                    OK [built from source and cached in 1.4s]
Successfully installed 1 package in 1.4 seconds.
- Project '~/R/gdf6y-analysis-scripts/rnaseq-Mmu' loaded. [renv 1.0.3]
The following package(s) have broken symlinks into the cache:
- annotate
- AnnotationDbi
...
...
...
- zip
- zlibbioc
Use `renv::repair()` to try and reinstall these packages.

- One or more packages recorded in the lockfile are not installed.
- Use `renv::status()` for more details.
```


```
> renv::repair()
# Library cache links --------------------------------------------------------
# Downloading packages -------------------------------------------------------
- Downloading annotate from Bioconductor ...    OK [file is up to date]
- Downloading Biobase from Bioconductor ...     OK [file is up to date]
...
...
...
- Downloading ProtGenerics from Bioconductor ... OK [8.5 Kb in 0.64s]
- Downloading GenomicFeatures from Bioconductor ... OK [1.3 Mb in 2.0s]
- Downloading openxlsx from CRAN ...            OK [file is up to date]
- Downloading zip from CRAN ...                 OK [file is up to date]
Successfully downloaded 150 packages in 440 seconds.

The following package(s) will be installed:
- annotate               [1.72.0]
- AnnotationDbi          [1.56.2]
- AnnotationFilter       [1.18.0]
- AnnotationHub          [3.2.2]
...
...
...
- zip                    [2.3.0]
- zlibbioc               [1.40.0]
These packages will be installed into "~/R/gdf6y-analysis-scripts/rnaseq-Mmu/renv/library/R-4.1/x86_64-pc-linux-gnu".

Do you want to proceed? [Y/n]: y
```


```
# Installing packages --------------------------------------------------------
- Installing BiocVersion ...                    OK [built from source and cached in 0.5s]
- Installing BiocGenerics ...                   OK [built from source and cached in 2.4s]
...
...
...
- Installing zip ...                            OK [built from source and cached in 5.7s]
- Installing openxlsx ...                       OK [built from source and cached in 29s]
Successfully installed 150 packages in 57 minutes.

# Package sources ------------------------------------------------------------
- All installed packages appear to be from a known source.
```
