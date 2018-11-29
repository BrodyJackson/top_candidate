# **Top Candidate Test**

[![Build Status](https://travis-ci.org/joemccann/dillinger.svg?branch=master)](https://travis-ci.org/joemccann/dillinger)

--- 


# Table of Contents
1. [Introduction](#introduction)
2. [Installation](#install)
3. [Function Documentation](#functiondocumentation)
   + [Experimental Data Import](#experimentaldataimport)
   + [ortholog_filtering](#orthologfiltering)
   + [Find_top_candidates](#findcandidates)
   + [compare_top_candidates](#compare)
4. [Demo usage](#demousage)
5. [Data format requirements](#dataformat)

--- 
 <a name="introduction"></a>
# Introduction

The top candidate test is a method created by Dr. Samuel Yeaman at the University of Calgary to isolate specific gene candidates from large sets of genomic data. Top candidates are selected by comparing the number of accurances to some measure of significance, then searching for those above an aribitrary outlier threshold. This tool was initially developed to find signatures of local adaptation, by finding outliers in the association between candidate density and phenotype-treatment corellation. In its original use the method identified putative locally adapted loci within species of Pine and Spruce trees. The correlation between various SNP's (phenotype) and environment (treatment) was compared to the density of each SNP, and outliers from this association were marked as top candidates. This package attempts to generalize the procedure, and supply useful methods for those wanting a similar way to parse candidates from their own experimental data. 

---
<a name="install"></a>
# Installation
The topCandidateTest package can be downloaded directly from github by using the following R commands

```
library(devtools)
install_github('topCandidateTest', 'github_BrodyJackson')
```

---

<a name="functiondocumentation"></a>
# Function Documentation 
<a name="experimentaldataimport"></a>
### experimental_data_import
###### Allows for input of experimental info about a species
  -  `Parameter`:  **overall_info** </br> 
     required table formatted data containing all the experimental info about a species. Defaults to NULL
   - `Parameter`: **get_only_annotated** </br>
     will get only the SNP's annotated to a gene if true, defaults to TRUE
  - `Parameter`: **annotatedColumn** </br>
     the column which you want to ensure that each value is annotated to. As an example this could be the contig column if you only wanted SNPS annotated to a gene. Defaults to NULL
  - `Keywords`: **species** 
  - `Examples`:
      ```
      experimental_data_import("all_spruce_data_pvalues.txt", TRUE, "tcontig")
      ```
<a name="orthologfiltering"></a>
### ortholog_filtering
###### Function to filter data for species entered
  -  `Description`: </br>
     will get the unique names of value asked for, such as genes,in species then find which ones have above a certain number of occurances
   - `Parameter`: **species_data** </br>
     data to evaluate. For example, this could be the contigs present in a species. Defaults to NULL
  - `Parameter`: **minimum** </br>
     the minimum number of occurences a value needs to have in this data to be included. For example, number of SNP's needed for a gene. Defaults to 3
  - `Parameter`: **columnName** </br>
     the column name which has the values to evaluate. For example, this could be the column of contigs. Defaults to NULL
  - `Keywords`: **ortholog** 
  - `Examples`: 
      ```
      ortholog_filtering(pine_tcontig, 3, "tcontig")
      ```
<a name="findcandidates"></a>
### find_top_candidates
###### Function to filter data for species entered
  -  `Description`: </br>
   - `Parameter`: **species** </br>
     data list containing overall data in first position and orthologs data to evaluate in second. These values are previously found using experimental_data_import method. Defaults to NULL
  - `Parameter`: **quantiles** </br>
     the quantiles being tested
  - `Parameter`: **binom** </br>
     the binomial cuts to test
  - `Parameter`: **candidateColumnName** </br>
     the name of the column containing all the candidates that will be tested. Defaults to NULL
  - `Parameter`: **testVariable** </br>
     the column name holding the results for a certain test. For example, this could be the results for a certain environmental location
  - `Parameter`: **plotFileName** </br>
     the name of the pdf file you want to plot in, if NULL then top candidates will not be plotted, defaults to null
  - `Return`: </br>
     returns a list containing data frames holding the top candidate results for each quantile
  - `Keywords`: **topCandidate** 
  - `Examples`: 
      ```
      find_top_candidates(pine, the_quantiles, binom_cuts, "tcontig", "MCMT", "pine_top_candidates.pdf")
      ```


<a name="compare"></a>
### compare_top_candidates
###### Function to filter data for species entered
  -  `Description`: </br>
   - `Parameter`: **species_one_candidates** </br>
     candidates list containing data frames holding top candidate values for various quantiles of species one. Defaults to NULL
  - `Parameter`: **species_two_candidates** </br>
     list containing data frames holding top candidate values for various quantiles of species two. Defaults to NULL
  - `Parameter`: **quantiles** </br>
     the quantiles being tested
  - `Parameter`: **binom** </br>
     the binomial cuts to test
  - `Parameter`: **orthologs** </br>
     data frame containing the one to one orthologs of the value being tested for the two species
  - `Parameter`: **ortho_in_both** </br>
     data frame containing the orthologs that are present in both the species being tests
  - `Parameter`: **species_one_name** </br>
     the name of species one (that is being used as column title)
   - `Parameter`: **species_two_name** </br>
     the name of species two (that is being used as column title)
   - `Parameter`: **plotFileName** </br>
     the name of the pdf file you want to plot in, if NULL then top candidates will not be plotted, defaults to null
  - `Return`: </br> 
     returns the probability of getting number_overlap or more top candidates
  - `Keywords`: **compare** 
  - `Examples`: 
    ```
    compare_top_candidatescompare_top_candidates(pine_outliers, spruce_outliers, the_quantiles, binom_cuts, oneToOne, ortho_in_both, "pine", "spruce", "top_candidates_compared.pdf")
    ```


---
<a name="demousage"></a>
# Demo Usage

The files present in the demo folder of the repository can be used to run testScript.R, which is what this demo will go over

The first step is to install the package for use in your script
```R
install("topCandidateTest")
library(topCandidateTest)
```
Next we want to read in the original data between pine and spruce orthologs, including both the one-to-one and multi-to one correspondences
```R
allOrthologs <- read.table ("all_orthologs_pine_spruce.txt",T)
oneToOne  <- read.table ("one_to_one_orthologs.txt", T)
```
Then we will use the experimental_data_import function on both species which returns a data frame that we can use to access both the overal experimental results from the data, as well as the results for only SNP's that are annoated to a gene
```R
pine <- experimental_data_import("all_pine_data_pvalues.txt", TRUE, "tcontig")
overall_pine <- data.frame(pine[[1]])
pine_tcontig <- data.frame(pine[[2]])

spruce <- experimental_data_import("all_spruce_data_pvalues.txt", TRUE, "tcontig")
overall_spruce <- data.frame(spruce[[1]])
spruce_tcontig <- data.frame(spruce[[2]])
```
We then filter the data we have about SNP's annotated to a gene, and get back only those that have more than 3 occurances
```R
filtered_orthologs_pine <- ortholog_filtering(pine_tcontig, 3, "tcontig")
filtered_orthologs_spruce <- ortholog_filtering(spruce_tcontig, 3, "tcontig")
```
Orthologs that are present in the data tables of both species must then be found
```R
int1 <- oneToOne[oneToOne$pine %in% filtered_orthologs_pine,]
ortho_in_both <- int1[int1$spruce %in% filtered_orthologs_spruce,]
```
We can now select the quantiles and binomial cuts that our candidates should be tested with.
```R
the_quantiles <- c (0.001,0.005,0.01,0.05)

binom_cuts <- c(0.99,0.999,0.9999,0.99999,0.999999,0.9999999,0.99999999,0.999999999)
```
The top candidates can now be found by using the formatted data, if a pdf file is supplied as a parameter then plots of top candidates can be found there
```R
pine_outliers <- find_top_candidates(pine, the_quantiles, binom_cuts, "tcontig", "MCMT", "pine_top_candidates.pdf")

spruce_outliers <- find_top_candidates(spruce, the_quantiles, binom_cuts, "tcontig", "MCMT", "spruce_top_candidates.pdf")
```
The top candidates from both species can now be compared, and the probability of getting a certain number of candidates that overlap in both species will be returned. If a pdf file is supplied then plots of these comparisons can be found there
```R
probability <- compare_top_candidates(pine_outliers, spruce_outliers, the_quantiles, binom_cuts, oneToOne, ortho_in_both, "pine", "spruce", "top_candidates_compared.pdf")

print(probability)
```
