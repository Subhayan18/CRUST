This package is created to offer a consolidated tool for clonal deconvolution while using bulk sequencing data when a tumor is multiregionally or 
multi-temporally sampled. The package addresses several shortcomings while subclonally reconstructing several samples together that show up quite often 
due to either sample quality, processing, sequencing or inferential artifacts.

## Installation
*CRUST* borrows clustering programs and other supporting facades from several other dependencies. Most of which has a continuing support from *CRAN*. Although more often than not one might find unable to import one or more dependencies due to lack of support from its maintainer. If such consern arises please refer to the package manual to 
find a list of direct and suggested dependencies which may help to resolve the issue.

```{r}
install.packages("remotes")
remotes::install_github("ShixiangWang/copynumber")
remotes::install_github("Subhayan18/CRUST")
```

If no error message has popped up, we are now ready to load the package in global environment.

```{r}
require(CRUST)
```

## Preliminary usage

*CRUST* comes with one simulated tetraploid tumor data that has eight different samples. This is mostly used for demonstration purpose. 
We will use this data to perform our first clonal deconvolution. A set of example user input is given below but there are lot of other choices 
in methods and inputs that can be provided.

The *cluster.doc* function needs to be called for the analysis which requires the user to declare the colum numbers of the sample names, variant allele frequencies 
and among other options, a method to select the optimum number of clusters and a method to perform clustering. The last two options have a default input if left 
undeclared.

As the variant allele frequency or VAF is the currency to any deconvolution program, it is paramount to ensure its integrity prior to analysis. A major concern 
is thus to adjust the VAFs according to the segmental copy numbers of each sample. To this end CRUST only analyzes segments that have same allelic makeup across 
samples. For example, in our imaginary tetraploid tumor the allelic makeup can assume the following states:

0+4 : 0 mutatant and 4 wildtypes

1+3 : 1 mutatant and 3 wildtypes

2+2 : 2 mutatant and 2 wildtypes

3+1 : 3 mutatant and 1 wildtypes

4+0 : 4 mutatant and 0 wildtypes

Instead of adjusting the VAFs corresponding to their respective segmental copy numbers *CRUST* analyzes each distinct allelic composition in separate runs. The user is prompted to declare this allelic make up in the beginning of each run. In addition to this *CRUST* also asks the user to make a visual inspection of the VAF dot plot and provide 
a intuitive solution to the problem. It does not uses this intuitive solution to infer cluster separation rather builds upon this and later in the analyses provides the user 
with options to modify solutions.

```{r}
## let's have a look at the data itself
?test.dat
head(test.dat)

res.1 <- cluster.doc(test.dat, sample = 1, vaf = 2, 
         optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 2
## Would you like to see my suggestion instead? (yes/no): no
```

**Figure 1** shows how the distribution of variant allele frequencies look for the 8 simulated samples. Depending on the structure of the spread it is concievable that the allelic segmentation is a balanced 2 + 2

**Figure 2** shows the changes in *Bayesian Information criteria* (BIC) estimated based on expectation and variance of the clustering fit.

**Figure 3** shows the clustered samples here depict the distribution of clonal and sub-clonal variants.

![](https://github.com/Subhayan18/CRUST/blob/gh-pages/test.dat.1.png)

![](https://github.com/Subhayan18/CRUST/blob/master/source/test.dat.2.png)

<p><img alt="" src="https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.3.png" /></p>

## Scaling

This data is rather clean and the clusters are more or less obviously identifiable with a visual inspection. But real data tend to be much more noisy and the clonal clusters  are unassumingmore often than not. Let's use the data set from the neuroblastoma patient sample *ES* to see how messy real data can be.

```{r, eval=FALSE, echo=TRUE}
es<-read.table("ES_all_2+2.txt",header=T,stringsAsFactors=F)

## Let's check the VAF distribution of the samples
ggplot(es, aes(x=sample, y=mut, col=as.factor(sample))) + geom_point()
```

<center>

![VAF distribution for ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.4.png)

</center>

What happend if we deconvolute this data assuming there are two clones and one subclone?

```{r, eval=FALSE, echo=TRUE}
res.2 <- cluster.doc(es, sample = 1, vaf = 2, 
                   optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## What is the suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 1
```

<center>

![Clonal deconvolution for ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.5.png)

</center>

As it is evident that the VAF distribution here is quite varied among samples, let's now use *Probabilistic Quotient Normalization* with the *Cancer Cell Fractions* (CCF)

```{r, eval=FALSE, echo=TRUE}
es.sc<-seqn.scale(es,2,3)

## Let's check the scaled VAF distribution of the samples

ggplot(es.sc, aes(x=sample, y=scaled.vaf, col=as.factor(sample))) + geom_point()
```

<center>

![VAF distribution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.6.png)

</center>

This step has pretty much rescaled the VAFs for the samples that had a relatively lower CCFs. Let's try the deconvolution now.

```{r, eval=FALSE, echo=TRUE}
res.2 <- cluster.doc(es.sc, sample = 1, vaf = 3, 
                   optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## What is the suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 2
## Would you like to see my suggestion instead? (yes/no): no
```

<center>

![Clonal deconvolution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.7.png)

</center>

Can you spot the difference in the subclonal distributions for some of the samples now?
*Hint*: observe how more data points become clonal after scaling for the first few samples.

Now let us assume we are not very satisfied as how things stand for **sample_6** and **sample_7**. Instead of **2 clonal and 1 sub-clonal clouds** we would rather see **2 clones and 2 subclones** for both samples.

Here's how we can make that happen:

```{r, eval=FALSE, echo=TRUE}
## Let's have a look at the function we can use for this:

?cluster doubt

## Now to the fitting:

res.3 <- cluster.doubt(res.2,1,3,c("sample_6","sample_7"),c(2,2,2,2))
```

<center>

![User rectified clonal deconvolution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.8.png)

</center>

## Estimation of allelic composition

When the allelic make up (Copynumber data from SNP array) is unavilable to the user, it can be estimated given the sequence reads from the constitutional DNA is also present. This can generally be obtained from a .vcf file before the variant calling.

```{r, eval=FALSE, echo=TRUE}
## A user provided .vcf file must contain data from one tumor sample
## and a corresponding normal tissue sample

m16 <- vcfR::read.vcfR("tumor_data.vcf")
CN.est <- AlleleComp(data=m16, AD = "AD", method = "apriori")
```

## Auxiliary functions

**Mutect2** is a popular variant caller that can be used to obtain quality controlled variant calls. Those calls can be tabulated with the **GATK** extension [**VariantsToTable**](https://gatk.broadinstitute.org/hc/en-us/articles/360042476292-VariantsToTable). **CloneStrat** offers an extension that can convert such an output to a data file compatible with the internal functions described.

```{r, eval=TRUE, echo=TRUE}
## Example whole exome seq data is provided in this repository
## Let's query the function to be used:

?mutect2.qc

## Time to see it in action

WES <- readxl::read_excel("WES.xlsx")
sample.name <- c("X1","X2","X3","X4","X5")
CS.dat <- mutect2.qc(WES,sample.name)
```

**A test** is provided to check the goodness of the cluster fit. This will indicate existence of outliers.
```{r, eval=FALSE, echo=TRUE}
CS.test<-T.goodness.test(es)$rej

## The variants that can only belong to one clone is shown
```
### Coming soon

Functions to plot copy number estimation and figure out allelic composition. I am also testing a new method for estimating copy numbers which gives user more freedom to tweak.

More features will be added gradually. If you'd like to see a specific feature incorporated in `CloneStrat`, send me 
a request [here](https://htmlpreview.github.io/?https://github.com/Subhayan18/CloneStrat/blob/master/footer.html)
