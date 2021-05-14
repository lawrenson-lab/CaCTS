CaCTS Tutorial
================

## R Markdown

**Ca**ncer **C**ore **T**ranscription factor **S**pecificity (CaCTS) is
an algorithm to identify candidate MTFs using pan-cancer RNA-sequencing
data from The Cancer Genome Atlas. In this page we describe the
procedures to perform CaCTS analysis. For more results details, please
check our [biorixv](https://www.biorxiv.org/content/10.1101/839142v1)
preprint .

# 1\. Installation

Install all packages in the latest version of
[R](https://www.r-project.org/).

### Option 1

Downaload
[CaCTS](https://github.com/lawrenson-lab/CaCTS_zip/raw/master/CaCTS_1.0.tar.gz)
package and run the following command in a *Terminal*

``` r
R CMD INSTALL CaCTS_1.0.tar.gz
```

### Option 2

First, install the devtools package.

``` r
install.packages("devtools")
```

Load the devtools package.

``` r
library(devtools)
```

Use install\_github as follows: Load the devtools package.

``` r
install_github("lawrenson-lab/CaCTS")
```

Load the CaCTS
    package.

``` r
library(CaCTS)
```

    ##  ===================================================================
    ##          ___        ___ _____  ___             _                    
    ##         |          |      |   |            _  | |                   
    ##         |    ___   |      |   |___        | | | | _                 
    ##         |     __|  |      |       |       | |_| || |                
    ##         |___ |__|  |___   |    ___|       |___   __|                
    ##                                               | |                   
    ##                                               | |                   
    ##  -------------------------------------------------------------------
    ##    Finding Cancer Core Transcription factor Specificity (CaCTS)     
    ##        Version:1.0
    ##  ===================================================================

## Pre-processing input files

Basically, CaCTS algorithm requires two main data structure inputs: the
*expression matrix* and the *annotation table*.

``` r
# Exoression data
load("data/TCGA.RNA.Rda")
TCGA.RNA[1:4, 1:4]
```

    ##   gene_symbol   gene_id TCGA-OR-A5J1-01A-11R-A29S-07
    ## 1           ? 100130426                       0.0000
    ## 2           ? 100133144                       3.2661
    ## 3           ? 100134869                       3.9385
    ## 4           ?     10357                     149.1350
    ##   TCGA-OR-A5J2-01A-11R-A29S-07
    ## 1                       0.0000
    ## 2                       2.6815
    ## 3                       8.9948
    ## 4                      81.0777

``` r
# Annotations
annot.table <- read.delim("data/SuppTable1-34-TCGAID.txt")
colnames(annot.table)[3] <- "sample.id"
colnames(annot.table)[2] <- "group.name"
head(annot.table)
```

    ##   Cancer group.name                    sample.id
    ## 1    ACC        ACC TCGA-OR-A5J1-01A-11R-A29S-07
    ## 2    ACC        ACC TCGA-OR-A5J2-01A-11R-A29S-07
    ## 3    ACC        ACC TCGA-OR-A5J3-01A-11R-A29S-07
    ## 4    ACC        ACC TCGA-OR-A5J5-01A-11R-A29S-07
    ## 5    ACC        ACC TCGA-OR-A5J6-01A-31R-A29S-07
    ## 6    ACC        ACC TCGA-OR-A5J7-01A-11R-A29S-07

``` r
annot.table$group.name <- paste0(annot.table$Cancer, "-", annot.table$group.name)
head(data.frame(table(annot.table$group.name)), 20)
```

    ##                                      Var1 Freq
    ## 1                                 ACC-ACC   78
    ## 2                               BLCA-BLCA  404
    ## 3                               BRCA-BRCA 1083
    ## 4                               CESC-CESC  301
    ## 5                               CHOL-CHOL   36
    ## 6                               COAD-COAD  441
    ## 7                               DLBC-DLBC   48
    ## 8      ESCA-Esophagus Adenocarcinoma  NOS   88
    ## 9  ESCA-Esophagus Squamous Cell Carcinoma   94
    ## 10                                GBM-GBM  156
    ## 11                              HNSC-HNSC  514
    ## 12                              KICH-KICH   65
    ## 13                              KIRC-KIRC  515
    ## 14                              KIRP-KIRP  285
    ## 15                              LAML-LAML  173
    ## 16                                LGG-LGG  514
    ## 17                              LIHC-LIHC  368
    ## 18                              LUAD-LUAD  510
    ## 19                              LUSC-LUSC  487
    ## 20                              MESO-MESO   87

In this step, from all Pancancer expressed genes, we are selecting only
TF genes. Published lists of transcription factors (TF) were retrieved
from Saint-André et al., 2016 and Lambert et al., 2018. Merging both
lists we created a catalogue of 1,671 unique TFs, of which 1,578 were
expressed in the pancan data set.

``` r
TF.list = read.delim("data/merged.list.1671.TFs.txt", sep = "\t")

f.TCGA.RNA = TCGA.RNA[which(TCGA.RNA$gene_symbol %in% as.character(TF.list$NameTF)),]
dim(f.TCGA.RNA)
```

    ## [1]  1578 11071

``` r
aux = which(!as.character(TF.list$NameTF) %in% f.TCGA.RNA$gene_symbol)
write.table(TF.list[aux,],file = "data/non-expressed-TFs.txt", quote = F, row.names = F)

rownames(f.TCGA.RNA) <- f.TCGA.RNA$gene_symbol
f.TCGA.RNA <- f.TCGA.RNA[,3:ncol(f.TCGA.RNA)]
f.TCGA.RNA[1:5, 1:4]
```

    ##       TCGA-OR-A5J1-01A-11R-A29S-07 TCGA-OR-A5J2-01A-11R-A29S-07
    ## ADNP2                      361.671                      693.732
    ## ADNP                      1452.450                     2917.060
    ## AEBP1                      471.182                     1758.690
    ## AEBP2                      842.459                     1689.440
    ## AFF3                      1782.900                     2962.550
    ##       TCGA-OR-A5J3-01A-11R-A29S-07 TCGA-OR-A5J5-01A-11R-A29S-07
    ## ADNP2                     350.1700                     329.9770
    ## ADNP                     1506.1500                    2087.5300
    ## AEBP1                     251.8150                     493.4160
    ## AEBP2                    1224.1100                    1317.5800
    ## AFF3                       43.2529                      54.9961

``` r
dim(f.TCGA.RNA)
```

    ## [1]  1578 11069

To calculate a Jensen-Shannon Divergence (JSD) score for each TF for
each tumor type, we adjusted the normalized expression values, to
rescale all values to \>0.

``` r
delta1 = max(f.TCGA.RNA, na.rm = T) - min(f.TCGA.RNA, na.rm = T)
delta2 = max(f.TCGA.RNA, na.rm = T) - 0

f.TCGA.RNA.rs = (f.TCGA.RNA - min(f.TCGA.RNA, na.rm = T)) * delta1 / delta2
f.TCGA.RNA.rs[1:5, 1:4]
```

    ##       TCGA-OR-A5J1-01A-11R-A29S-07 TCGA-OR-A5J2-01A-11R-A29S-07
    ## ADNP2                     362.5552                     694.6164
    ## ADNP                     1453.3349                    2917.9460
    ## AEBP1                     472.0662                    1759.5752
    ## AEBP2                     843.3435                    1690.3251
    ## AFF3                     1783.7852                    2963.4360
    ##       TCGA-OR-A5J3-01A-11R-A29S-07 TCGA-OR-A5J5-01A-11R-A29S-07
    ## ADNP2                    351.05415                    330.86113
    ## ADNP                    1507.03498                   2088.41540
    ## AEBP1                    252.69908                    494.30025
    ## AEBP2                   1224.99478                   1318.46485
    ## AFF3                      44.13683                     55.88004

``` r
dim(f.TCGA.RNA.rs)
```

    ## [1]  1578 11069

## Step 1 - Representative sample

Next, for each tumor type, a representative sample was generated by
calculating the mean expression value for each
TF.

``` r
matrix.rep <- prepare_representaive_samples(expr.matrix = f.TCGA.RNA.rs, sample.descr = annot.table, save.file = F)
matrix.rep[1:4, 1:4]
```

    ##         ACC-ACC BLCA-BLCA  BRCA-BRCA CESC-CESC
    ## ADNP2  467.2289  649.4564   765.7106  700.8787
    ## ADNP  1898.5172 2455.4685  3606.1140 2393.4897
    ## AEBP1 1011.8737 6894.9661 15561.8556 4121.5774
    ## AEBP2 1245.3867  845.4508   632.1004  954.5149

## Step 2 - Calculatiog the JSD score

``` r
cancer = "OV-OV"
res.CaCTS <- run_CaCTS_score(matrix.rep, cancer)
```

    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |                                                                 |   1%
      |                                                                       
      |=                                                                |   1%
      |                                                                       
      |=                                                                |   2%
      |                                                                       
      |==                                                               |   2%
      |                                                                       
      |==                                                               |   3%
      |                                                                       
      |==                                                               |   4%
      |                                                                       
      |===                                                              |   4%
      |                                                                       
      |===                                                              |   5%
      |                                                                       
      |====                                                             |   5%
      |                                                                       
      |====                                                             |   6%
      |                                                                       
      |====                                                             |   7%
      |                                                                       
      |=====                                                            |   7%
      |                                                                       
      |=====                                                            |   8%
      |                                                                       
      |======                                                           |   8%
      |                                                                       
      |======                                                           |   9%
      |                                                                       
      |======                                                           |  10%
      |                                                                       
      |=======                                                          |  10%
      |                                                                       
      |=======                                                          |  11%
      |                                                                       
      |=======                                                          |  12%
      |                                                                       
      |========                                                         |  12%
      |                                                                       
      |========                                                         |  13%
      |                                                                       
      |=========                                                        |  13%
      |                                                                       
      |=========                                                        |  14%
      |                                                                       
      |=========                                                        |  15%
      |                                                                       
      |==========                                                       |  15%
      |                                                                       
      |==========                                                       |  16%
      |                                                                       
      |===========                                                      |  16%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |===========                                                      |  18%
      |                                                                       
      |============                                                     |  18%
      |                                                                       
      |============                                                     |  19%
      |                                                                       
      |=============                                                    |  19%
      |                                                                       
      |=============                                                    |  20%
      |                                                                       
      |=============                                                    |  21%
      |                                                                       
      |==============                                                   |  21%
      |                                                                       
      |==============                                                   |  22%
      |                                                                       
      |===============                                                  |  22%
      |                                                                       
      |===============                                                  |  23%
      |                                                                       
      |===============                                                  |  24%
      |                                                                       
      |================                                                 |  24%
      |                                                                       
      |================                                                 |  25%
      |                                                                       
      |=================                                                |  25%
      |                                                                       
      |=================                                                |  26%
      |                                                                       
      |=================                                                |  27%
      |                                                                       
      |==================                                               |  27%
      |                                                                       
      |==================                                               |  28%
      |                                                                       
      |===================                                              |  29%
      |                                                                       
      |===================                                              |  30%
      |                                                                       
      |====================                                             |  30%
      |                                                                       
      |====================                                             |  31%
      |                                                                       
      |=====================                                            |  32%
      |                                                                       
      |=====================                                            |  33%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |======================                                           |  34%
      |                                                                       
      |======================                                           |  35%
      |                                                                       
      |=======================                                          |  35%
      |                                                                       
      |=======================                                          |  36%
      |                                                                       
      |========================                                         |  36%
      |                                                                       
      |========================                                         |  37%
      |                                                                       
      |========================                                         |  38%
      |                                                                       
      |=========================                                        |  38%
      |                                                                       
      |=========================                                        |  39%
      |                                                                       
      |==========================                                       |  39%
      |                                                                       
      |==========================                                       |  40%
      |                                                                       
      |==========================                                       |  41%
      |                                                                       
      |===========================                                      |  41%
      |                                                                       
      |===========================                                      |  42%
      |                                                                       
      |============================                                     |  42%
      |                                                                       
      |============================                                     |  43%
      |                                                                       
      |============================                                     |  44%
      |                                                                       
      |=============================                                    |  44%
      |                                                                       
      |=============================                                    |  45%
      |                                                                       
      |==============================                                   |  45%
      |                                                                       
      |==============================                                   |  46%
      |                                                                       
      |==============================                                   |  47%
      |                                                                       
      |===============================                                  |  47%
      |                                                                       
      |===============================                                  |  48%
      |                                                                       
      |================================                                 |  48%
      |                                                                       
      |================================                                 |  49%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |=================================                                |  50%
      |                                                                       
      |=================================                                |  51%
      |                                                                       
      |=================================                                |  52%
      |                                                                       
      |==================================                               |  52%
      |                                                                       
      |==================================                               |  53%
      |                                                                       
      |===================================                              |  53%
      |                                                                       
      |===================================                              |  54%
      |                                                                       
      |===================================                              |  55%
      |                                                                       
      |====================================                             |  55%
      |                                                                       
      |====================================                             |  56%
      |                                                                       
      |=====================================                            |  56%
      |                                                                       
      |=====================================                            |  57%
      |                                                                       
      |=====================================                            |  58%
      |                                                                       
      |======================================                           |  58%
      |                                                                       
      |======================================                           |  59%
      |                                                                       
      |=======================================                          |  59%
      |                                                                       
      |=======================================                          |  60%
      |                                                                       
      |=======================================                          |  61%
      |                                                                       
      |========================================                         |  61%
      |                                                                       
      |========================================                         |  62%
      |                                                                       
      |=========================================                        |  62%
      |                                                                       
      |=========================================                        |  63%
      |                                                                       
      |=========================================                        |  64%
      |                                                                       
      |==========================================                       |  64%
      |                                                                       
      |==========================================                       |  65%
      |                                                                       
      |===========================================                      |  65%
      |                                                                       
      |===========================================                      |  66%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |============================================                     |  67%
      |                                                                       
      |============================================                     |  68%
      |                                                                       
      |=============================================                    |  69%
      |                                                                       
      |=============================================                    |  70%
      |                                                                       
      |==============================================                   |  70%
      |                                                                       
      |==============================================                   |  71%
      |                                                                       
      |===============================================                  |  72%
      |                                                                       
      |===============================================                  |  73%
      |                                                                       
      |================================================                 |  73%
      |                                                                       
      |================================================                 |  74%
      |                                                                       
      |================================================                 |  75%
      |                                                                       
      |=================================================                |  75%
      |                                                                       
      |=================================================                |  76%
      |                                                                       
      |==================================================               |  76%
      |                                                                       
      |==================================================               |  77%
      |                                                                       
      |==================================================               |  78%
      |                                                                       
      |===================================================              |  78%
      |                                                                       
      |===================================================              |  79%
      |                                                                       
      |====================================================             |  79%
      |                                                                       
      |====================================================             |  80%
      |                                                                       
      |====================================================             |  81%
      |                                                                       
      |=====================================================            |  81%
      |                                                                       
      |=====================================================            |  82%
      |                                                                       
      |======================================================           |  82%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |======================================================           |  84%
      |                                                                       
      |=======================================================          |  84%
      |                                                                       
      |=======================================================          |  85%
      |                                                                       
      |========================================================         |  85%
      |                                                                       
      |========================================================         |  86%
      |                                                                       
      |========================================================         |  87%
      |                                                                       
      |=========================================================        |  87%
      |                                                                       
      |=========================================================        |  88%
      |                                                                       
      |==========================================================       |  88%
      |                                                                       
      |==========================================================       |  89%
      |                                                                       
      |==========================================================       |  90%
      |                                                                       
      |===========================================================      |  90%
      |                                                                       
      |===========================================================      |  91%
      |                                                                       
      |===========================================================      |  92%
      |                                                                       
      |============================================================     |  92%
      |                                                                       
      |============================================================     |  93%
      |                                                                       
      |=============================================================    |  93%
      |                                                                       
      |=============================================================    |  94%
      |                                                                       
      |=============================================================    |  95%
      |                                                                       
      |==============================================================   |  95%
      |                                                                       
      |==============================================================   |  96%
      |                                                                       
      |===============================================================  |  96%
      |                                                                       
      |===============================================================  |  97%
      |                                                                       
      |===============================================================  |  98%
      |                                                                       
      |================================================================ |  98%
      |                                                                       
      |================================================================ |  99%
      |                                                                       
      |=================================================================|  99%
      |                                                                       
      |=================================================================| 100%

``` r
tail(res.CaCTS)
```

    ##          Name     value  LogValue
    ## res.821 SOX17 0.3954705 0.4028859
    ## res.395 HOXD3 0.3722604 0.4291532
    ## res.188  EMX2 0.3634673 0.4395346
    ## res.119 CTCFL 0.3597423 0.4440085
    ## res.472  LHX1 0.3446310 0.4626457
    ## res.985   WT1 0.1666036 0.7783157

``` r
p = visualize_scores(tf.scores = res.CaCTS, topn = 0, ncol = 1, rep.matrix = matrix.rep, query.name = cancer, filename = paste0(cancer, "-CaCTS-scores.pdf"), w=5, h=5)
p
```

![](README_files/figure-gfm/visualize-1.png)<!-- -->

The final candidate MTF list for each cancer type was defined by
considering the intersection of the 5% most highly expressed TFs in each
tumor type, and the top 5% TFs when ranked by JSD
score.

``` r
filtered <- filter_by_expression_rank(rep.matrix = matrix.rep, tf.scores = res.CaCTS, query.name = cancer, pn=0.05, pnE=0.05)
```

    ## Filtering TFs by expression rank...

    ## File saved as: /home/abraaodem/Code/MTF/CaCTS/CaCTs_res/OV-OV-CaCTS-scores-expression-filtered.txt

``` r
filtered
```

    ##       Name     value  LogValue JSD.rank Expr.mean Expr.rank
    ## 12     WT1 0.1666036 0.7783157        1  5364.437        17
    ## 2     EMX2 0.3634673 0.4395346        4  2874.900        55
    ## 10   SOX17 0.3954705 0.4028859        6  3036.156        46
    ## 5    MEIS1 0.4361167 0.3603972        9  4125.429        25
    ## 1  BHLHE41 0.4867695 0.3126766       16  5811.158        13
    ## 7     PAX8 0.4944187 0.3059051       17  6817.798         9
    ## 3     ESR1 0.5004109 0.3006732       18  2811.067        60
    ## 14  ZNF503 0.5198643 0.2841100       28  2411.554        74
    ## 4    MECOM 0.5394078 0.2680828       35  3053.090        45
    ## 11   TGIF2 0.5543643 0.2562047       46  2305.685        77
    ## 6    NR2F6 0.5624774 0.2498949       58  3193.932        41
    ## 8     PBX1 0.5696953 0.2443573       64  4244.521        24
    ## 9   PLSCR1 0.5753966 0.2400327       76  2843.742        58
    ## 13  ZNF217 0.5759398 0.2396229       78  3122.163        43

To predict the candidate MTFs for the other cancer types, the same
procedures on **Step 2** can be performed by changin the cancer of
interest.

## References

D’Alessio, A.C., Fan, Z.P., Wert, K.J., Baranov, P., Cohen, M.A., Saini,
J.S., Cohick, E., Charniga, C., Dadon, D., Hannett, N.M., et al. (2015).
A Systematic Approach to Identify Candidate Transcription Factors that
Control Cell Identity. Stem Cell Reports 5, 763–775.

Lambert, S.A., Jolma, A., Campitelli, L.F., Das, P.K., Yin, Y., Albu,
M., Chen, X., Taipale, J., Hughes, T.R., and Weirauch, M.T. (2018). The
human transcription factors. Cell 175, 598–599.

Saint-André, V., Federation, A.J., Lin, C.Y., Abraham, B.J., Reddy, J.,
Lee, T.I., Bradner, J.E., and Young, R.A. (2016). Models of human core
transcriptional regulatory circuitries. Genome Res. 26, 385–396.
