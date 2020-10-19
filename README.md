# The *FBCwPlaid* user's guide #
2020/10/19 9:11:34   

## Function description ##
***FBCwPlaid*: A Functional Bi-clustering Analysis of Epi-transcriptome Profiling Data via a Weighted Plaid Model**

*FBCwPlaid*  based on Lagrange multiplier method  to discover the potential functional patterns. The model seeks for one bi-cluster each time. Thus, the goal of each time turns into a binary classification problem. It initializes model parameters by *k*-means clustering, and then updates the parameters of the Plaid model. To address the issue that site expression level determines methylation level confidence, it uses RNA expression levels of each site as weights to make lower expressed sites less confident. *FBCwPlaid* also allows overlapping bi-clusters, indicating some sites may participate in multiple biological functions.


## Preparation before using *FBCwPlaid*
### Data preparation ###
As is known, MeRIP-seq data profiles the m<sup>6</sup>A epi-transcriptome by **input data** and **IP data**.Therefore, before using FBCwPlaid, it is required to prepare the required IP samples as well as input samples. According to the sample data of **IP** and **input**, the **methylation level** and **expression level** of m<sup>6</sup>A sites were calculated.

Be careful, the two sets of data fed into *FBCwPlaid* should be **numerical matrices (or dataframe)**.


## Instructions for using *FBCwPlaid* ##
### Run *FBCwPlaid* to get the concrete potential functional patterns. ###

The operation of plaid should be set with the following input parameters.

**FPKM.IP:** IP sample data.

**FPKM.input:** input sample data.

**Methylation.level:** RNA methylation level.

**Expression.level:** RNA expression level.

**optimization:** Logical variables. If TRUE, the enrichment constraint module is enabled.

**GENES.CORRES.SITES,:** The gene corresponding to each m<sup>6</sup>A site.

**GENE.ID.TYPES:** The type of gene ID in GENES.CORRES.SITES. Four types are supported: "ENTREZID", "SYMBOL", "ENSEMBL" and "GENENAME".

**exponent.num:** The value of exponential power obtained by enrichment constraint.

**kmeans.startup:** The number of iterations of *k*-means.

**max.layers:** The maximum number of patterns allowed to be found.

**iter.layer:** The number of iterations to find each pattern.

**iter.bin:** The number of iterative searches in binary.

**backfitting.num:** The number of times of back fitting.

**verbose:** Logical variables. If TRUE, additional information about progress is printed.

***Note:** In general, it is unknown how many patterns there are in the data, so **max.layers** suggests setting it to a large number. Eventually, FBCwPlaid can automatically determine how many patterns are based on the decision conditions.*

Based on the above input, you can run *FBCwPlaid* with the following code:

#### 1. Open the enrichment constraint module of *FBCwPlaid*:

        exponent <- FBCwPlaid(Methylation.level = data, Expression.level = weight, max.layers = 10, 
    		            optimization = TRUE, GENES.CORRES.SITES = gene_id, GENE.ID.TYPES = "ENTREZID", 
    		            kmeans.startup = 3, iter.layer = 20, iter.bin = 10, backfitting.num = 3, verbose = "FALSE")
        # Output of FBCwPlaid parameter optimization result:
	  # exponent: The best choice of exponential power. 

	
#### 2. Running *FBCwPlaid* under the best exponential power obtained, a series of patterns are obtained:

        bicluster <- FBCwPlaid(Methylation.level = data, Expression.level = weight, max.layers = 10, 
    		               optimization = FALSE, exponent.num = exponent, kmeans.startup = 3, 
			       iter.layer = 20, iter.bin = 10, backfitting.num = 3, verbose = "FALSE")
        # Output items for running FBCwPlaid:
	# bicluster: The final biclustering result. bicluster has a biclust class that can be called by other functions, such as visualization using the *drawHeatmap* function in the R biclust package.
	# mu.all: The background of the overall data.
	# mu.rec: The background value of each pattern.
	# alpha.rec: The alpha value corresponding to the sites in each pattern.
	# beta.rec: The beta value corresponding to the conditions in each pattern.
	# CPS.rec: Significant CPS value of each pattern.



## Contact ##
Please contact the maintainer of *FBCwPlaid* if you have encountered any problems:

**Shutao Chen:** shutao.chen@cumt.edu.cn
