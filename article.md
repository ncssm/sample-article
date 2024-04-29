---
title: "An Analysis on Genomic Correlation for Gallstone Susceptibility"
abstract: |
  Gallstones, typically benign and harmless hardened deposits of digestive fluids, present in the gallbladder can cause severe painful complications if left untreated and can lead to removal surgery. (@noauthor_gallstones_2021) Quantitative Trait Loci (QTL) analysis can be used to find potential genetic links through the analysis of logarithm of the odds (LOD) scores which can indicate a possible connection between those loci on the mouse chromosome and phenotypic presentation of a trait (@myles2008quantitative) linked to gallstone susceptibility such as weight, presence, severity and liver weight. Analyses were performed using the R/QTL package on a cohort of mice, in an intercross breed and fed a high-fat diet. (@lyons1) Analyses were compared between augmented data sets (to possibly prevent overfitting) and an analysis run on a non-augmented data set.
---

## Introduction
The gallbladder is an organ in the upper right portion of the abdomen, directly below the liver, that releases bile, a fluid that the liver produces, that digests fats. Bile is a solution of cholesterol, bilirubin, and bile salts. (@clevelandclinic_gallbladder:_nodate) There are two primary categories of gallstones, cholesterol and pigment stones. Cholesterol stones form due to a lack of balance of cholesterol, bilirubin, and bile salts in the bile. It can form due to excess bilirubin or cholesterol or a lack of bile salts. The cause of pigment stones is currently unknown, but they tend to develop in patients already suffering from cirrhosis, biliary tract infections, and hereditary blood disorders such as sickle cell anemia. (@noauthor_gallstones_2021)

While some traits are very clearly and individually linked to a single particular spot on the genome, most traits are inherently complex and thus, there are multiple locations on the genome that can influence the manifestation of that trait. Gallstone disease and susceptibility fall under that umbrella of traits. Gallstone susceptibility is a multifaceted trait influenced by both genetic and environmental factors and represents a significant health concern with considerable variability in its occurrence among individuals. The development of gallstones is associated with a complex interplay of genetic predisposition and lifestyle factors such as diet and obesity. (@noauthor_gallstones_2021)

Quantitative Trait Loci (QTL) analyses can be used to find multiple locations on the genome with high logarithm of the odds (LOD) scores which can indicate a possible correlation or causation between the two. QTL analysis uses statistical methods to link quantitative phenotypic traits to genetic markers on the chromosome to try to genetically explain extremely complex phenotypes. (@myles2008quantitative).

We seek to uncover QTLs that may harbor candidate genes influencing gallstone formation.

The findings from this data-driven approach not only contribute to our understanding of the genetic determinants of gallstone susceptibility but also pave the way for potential insights into personalized preventive strategies and therapeutic interventions.

## Data set(s)

### Setup

In this study, we used the R/qtl package in R (@broman2003r), developed by Bromen et. al as well as Ritsert Jansen’s MQM method (@arends2010r) to perform our data analysis. R is a programming language and open-source software environment specifically designed for statistical computing and data analysis. Widely used by statisticians, data scientists, and researchers, R provides a comprehensive suite of tools for data manipulation, statistical modeling, visualization, and the development of custom analytical workflows. R/qtl is an R package designed for conducting Quantitative Trait Locus (QTL) analysis, a statistical method used to identify genetic loci associated with variations in quantitative traits. Developed for genetic mapping studies, R/qtl provides a range of tools for analyzing experimental cross populations, facilitating the detection and characterization of genomic regions influencing complex traits.


### Data collection

Data for this project was obtained from the Mouse Phenome Database at the Jackson Laboratory, a grant-funded resource that provides integrated genomic and phenomic data on behavioral, morphological, and physiological characteristics in mice. (@10.1007/s00335-023-10014-3) The Jackson Lab is an independent non-profit biomedical research lab that primarily conducts genomic research with mice.

```{figure} images/f2_breeding_chart.svg
:name: breeding_chart
:width: 500px
:align: center

Breeding Chart for Intercross Strains
```

@breeding_chart shows an in-depth chart, demonstrating visually how the mice are bred for effective analysis. This specific dataset, the Lyons1 data set, looks at plasma lipids and gallstone susceptibility in the F2 progeny of a DBA/2J x CAST/EiJ intercross. (@lyons1) There are two primary crosses of mice used for QTL analyses, intercross and backcross. As depicted by @breeding_chart, an intercross is characterized by two homozygous mice in the parental generation bred to produce heterozygous F1 or first-filial generation children. These F1 mice are then bred together to produce the F2 generation who are then utilized in experiments and studies. An interesting side note is that all of the AA, or homozygous dominant mice in the parental generation are deeply inbred and thus genetically identical, and so are all of the BB, or homozygous recessive mice in the parental generation. (@silver1995mouse [Chapter 3, Section 2])

Only male mice were utilized in the study. The animals had unrestricted access to both food and water and were housed in a temperature-regulated environment (71.6°F - 73.4°F approximately) which had a 14 hours of light and 10 hours of dark cycle. The animals were initially fed a low-cholesterol diet until the age of 6-8 weeks when they were switched to a lithogenic, high-cholesterol diet. This diet was composed of 15 percent butterfat, 1 percent cholesterol, 0.5 percent cholic acid, 2 percent corn oil, 50 percent sucrose, and 20 percent casein. All experimental protocols were approved by the Institutional Animal Care and Use Committees of The Jackson Laboratory and Harvard University. (@lyons2003quantitative)

### Data structure

```{figure} #genetic_map
:name: fig-genetic_map
:width: 500px
:align: center

Genetic Map for Markers on Mouse Chromosome
```

The data used was an F2 intercross with 278 individuals. There were 15 phenotypes and all 15 phenotypes had over 96.8 percent of the individual mice phenotyped. Mice have 20 chromosomes, 19 autosomes, and one sex chromosome, the X chromosome. There are 109 molecular markers in this data and the genetic map is shown above in @fig-genetic_map.

There was a 97.4 percent rate of genotyping, meaning this data is extremely complete and this is shown in @fig-missing_genotypes below. We can see that there is very little missing data which is 2.6 percent of total data missing according to the summary function in R/qtl.

Of the 15 phenotypes presented in the dataset, we chose to focus on four: a score on the severity of the gallstone, the number of gallstones measured, the weight of the gallstones, and the aggregates of the severity of cholesterol monohydrate crystals, which is a key indicator of gallstone development.

```{figure} #missing_genotypes
:name: fig-missing_genotypes
:width: 500px
:align: center

Missing Data Map
```

## Data Preparation and Modeling

### Data Preparation

Prior to running the analyses, we had to prepare the data further. We first completed a pairwise recombination factor plot to take a look at the physical distances between markers on the chromosome and ensure that they are accurate. We first estimate recombination fractions between markers within a genetic cross with the est.rf function. Recombination fractions are crucial in genetic mapping as they indicate the likelihood of genetic crossovers occurring between markers during the formation of gametes. These fractions are fundamental for constructing genetic maps, elucidating the distances between genetic markers, and identifying regions of the genome associated with specific traits through QTL analyses. We then generate a visual representation of the estimated recombination fractions which is essential in understanding the genetic linkage and physical distances between markers along the chromosomes, providing researchers with insights into the genetic architecture of traits and facilitating the identification of potential genomic regions influencing complex phenotypes. We then plotted our pairwise recombination scores and LOD scores in @fig-pairwise_recombination.

```{figure} #pairwise_recombination
:name: fig-pairwise_recombination
:width: 500px
:align: center

Pairwise Recombination Scores and LOD Scores
```

As evidenced by the lack of large red spots and a clean line, this data is clean and alright to use for further analysis.

### Exploratory Data Analysis

After cleaning our data set and ensuring quality control, we then began exploratory data analysis. The first step was using the R/ggplot2 library to explore trends within the phenotypic data. We used a correlation heat map to identify correlations between phenotypes and the manifestation of gallstones. Correlation heat maps are a visual representation of the coefficient of determination between various factors, or the r-squared value in a color-coded matrix. A value with an absolute value of 1 has a very strong correlation, either positive or negative. A value closer to 0 has less correlation and is more random. A positive number correlates to a positive correlation, meaning as the x-value increases, so does the y-value. A negative number correlates to a negative correlation, meaning as the x-value increases, the y-value decreases. This is shown below in @correlogram.

```{figure} images/image05.png
:name: correlogram
:width: 500px
:align: center

Correlation Heat Map
```

In this figure, we can see that there is a high positive correlation between multiple factors, particularly between the binary classification of the solidity of gallstones and the number of gallstones, the weight of gallstones and the number of gallstones, the weight of the gallstones and the binary solidity classification, and the presence of gallstones and the severity score.

We also ran a dendrogram heat map to analyze correlations as well as find hierarchical correlations between our phenotypic factors, shown in @fig-dendrogram. This map shows both a heat map to show correlations between various factors, similar to our previous graph, but also shows hierarchical relationships between our phenotypic factors. This helps us understand the degree of the relationship between the factors.

```{figure} #dendrogram
:name: fig-dendrogram
:width: 500px
:align: center

Dendrogram Hierarchical Heat Map
```

###  QTL Analyses

Following this, we can begin the setup for the QTL analysis. We first calculated the genotypic probabilities for individuals in the genetic cross by determining the likelihood of different genetic marker configurations based on specified parameters such as recombination step size, genotyping error probability, and the Haldane map function. We then simulated genotypic data for the markers in the genetic cross, incorporating factors like recombination, genotyping errors, and mapping functions. This is the primary first step we must do prior to running any analyses and determining the locus of interest.

Following this, we completed the same steps for each of our four factors to generate a main scan analysis as well as effect plots for the highest probability locus of interest as determined by LOD scores. We first used the function scanone with a normal model and the “em” method, which is the Expectation-Maximization method which estimates missing genotype probabilities in the genetic mapping analysis. While this is specifically excellent for data sets with missing phenotypic information, we find that it is still a robust method of analysis.

For each method, we then completed a permutation on the scan with 100 permutations to assess the significance of LOD peaks for each phenotype. We then assigned threshold values based on our permutation results with confidence intervals of 95 percent, 90 percent, and 63 percent respectively. After this, we plotted the results onto a main scan plot with colored lines representing our thresholds. We also ran a summary of the scan per phenotype and identified the most probable locus of interest per scan. We then used this location to identify a molecular marker in our data set and run an effect plot. These plots are particularly useful for understanding how genetic variation at specific loci influences the phenotypic variation in a quantitative trait. We can see how homozygous recessive or dominant, or heterozygous affects the manifestation of different phenotypic traits.

Completing this, we decided to explore augmented data to see if results run on augmented data on a total QTL analysis for the data set. To do this, we first created an augmented data set derived from our cross with a minprob of 0.1. This establishes a minimum probability threshold for considering the effects of additional markers or QTLs. This threshold, which influences the augmentation step, allows us to filter out less statistically significant QTL effects, refining the model and focusing on those with higher confidence. The choice of the minprob threshold serves as a key determinant in balancing sensitivity and specificity in the identification of quantitative trait loci, tailoring the analysis to the desired level of statistical rigor. We chose 0.1 because we wanted more statistically significant results rather than a broader overview with less statistically significant results. We then ran a geno.image on both the augmented cross and the original cross and compared the plots. Following this we took a scan of each cross, using the mqm scan for the augmented set and scanone for the original set, and found the peaks on each plot. We then found a molecular marker corresponding to each peak and compared it to each other. Following this, we used that marker we identified earlier as a cofactors and took another mqm scan with the cofactor of D18Mit64, the marker we identified. We then proceeded to plot all three main scans together on the same plot and compared the peaks.

## Results

### Expectation-Maximization Model
Using the Expectation-Maximization model, we generate 4 separate main scans with threshold lines at 95 percent confidence, 90 percent confidence, and 63 percent respectively.

#### Gall Count
As shown below in @fig-main_scan, for a count of gallstones present, we found two loci with a peak over 95 percent confidence.

```{figure} #main_scan
:name: fig-main_scan
:width: 500px
:align: center

Main Scan for Gall Count
```

There was 98 percent confidence in the correlation between c6.loc6 and phenotypic manifestation and 96 percent confidence in the correlation between c8.loc58 and phenotypic manifestation. We then found the correlated molecular markers for those two spots which were D6Mit46 and D8Mit88 respectively. Using those two markers, we then plotted an effect plot as shown below in @fig-effect_plot. As we can see in both figures, a homozygous DD genotype at both of these locations can correlate to gallstone susceptibility and a higher amount of gallstones while heterozygous DC and homozygous CC both present lower amounts of gallstones.

```{figure} #effect_plot
:name: fig-effect_plot
:width: 500px
:align: center

Effect Plot showing Allele vs Phenotypic Presentation
```

#### Gall Score

As shown below in @fig-main_scan_gall_score, for a score on how severe the gallstones were, we found one locus with a peak over or equal to 95 percent confidence.

```{figure} #main_scan_gall_score
:name: fig-main_scan_gall_score
:width: 500px
:align: center

Main Scan for Gall Score
```

There was 95 percent confidence of correlation between c2.loc52 and phenotypic manifestation. We then found the correlated molecular markers for this spot which was D2Mit94. Using this marker, we then plotted an effect plot as shown below in @fig-effect_plot_d2mit94. As we can see, a homozygous CC genotype at this location can correlate to a lower gallstone severity score while heterozygous DC and homozygous D both present higher scores. It can be said then that high gallstone severity is a dominant trait at this location in the genome.

```{figure} #effect_plot_d2mit94
:name: fig-effect_plot_d2mit94
:width: 500px
:align: center

Effect Plot showing Allele vs Phenotypic Presentation
```

#### Gall Weight
As shown below in @fig-main_scan_gall_weight, for a score on how severe the gallstones were, we found one locus with a peak over or equal to 95 percent confidence.

```{figure} #main_scan_gall_weight
:name: fig-main_scan_gall_weight
:width: 500px
:align: center

Main Scan for Gall Weight
```

There was 95 percent confidence in the correlation between c2.loc52 and phenotypic manifestation. We then found the correlated molecular markers for this spot which was D8Mit88. Using this marker, we then plotted an effect plot as shown below in @fig-effect_plot_d8mitt88. As we can see, a homozygous DD genotype at this location can correlate to a higher gallstone weight score while heterozygous DC and homozygous C both present lower weights. It can be said then that high gallstone weight is a recessive trait at this location in the genome.

```{figure} #effect_plot_d8mitt88
:name: fig-effect_plot_d8mitt88
:width: 500px
:align: center

Effect Plot showing Allele vs Phenotypic Presentation
```

#### Cholesterol Monohydrate Crystals, aggregates

As shown below in @fig-cholesterol_monohydrate_crystals, for a score on the cholesterol monohydrate crystals, we found one locus with a peak over or equal to 95 percent confidence.

```{figure} #cholesterol_monohydrate_crystals
:name: fig-cholesterol_monohydrate_crystals
:width: 500px
:align: center

Main Scan for Gall Weight
```

There was 95 percent confidence of correlation between c6.loc56 and phenotypic manifestation. We then found the correlated molecular markers for this spot which was D6Mit62. Using this marker, we then plotted an effect plot as shown below in @fig-effect_plot_d6mit62. As we can see, a homozygous DD genotype at this location can correlate to a lower cholesterol monohydrate crystal score while heterozygous DC and homozygous D both present higher scores.

```{figure} #effect_plot_d6mit62
:name: fig-effect_plot_d6mit62
:width: 500px
:align: center

Effect Plot showing Allele vs Phenotypic Presentation
```

### Augmented QTL Analyses Comparison

After generating an augmented dataset, we then used geno.image to plot both crosses respectively as shown below in @fig-plot_grid_og_ag.

```{figure} #plot_grid_og_ag
:name: fig-plot_grid_og_ag
:width: 500px
:align: center

Plot Grid of Original Genotype Data (left) and Augmented Genotype Data (right)
```

The genotypes CC, DC, and DD are displayed in the colors red, blue, and green, respectively. The white spaces represent missing data. As we can see, the augmented data is filled in much better, and there is no missing data. While running the summary function on the data, we see that nothing has changed in the augmented as compared to this original other than there being much more individuals in the data set (1343 as compared to 278). Additionally, the percent phenotyped remains approximately the same and the percent genotyped jumps up to 100 percent.

Next, we complete an mqmscan on the augmented dataset and a scanone on the original dataset and take a look at the maximum point on both of these. The augmented dataset has a peak at c18.loc5 which is 5 centimorgans on chromosome 18. The original dataset has a peak at c18.loc4 which is 8.46 centimorgans on chromosome 18. We then extract a marker for both of these positions which both are D18Mit64. We can then set that marker as a covariate and analyze for a new peak, with this marker as an additional variable.

We then plotted all three of these plots on the same map with green representing the original data, red representing the augmented data, and blue representing augmented data with the peak as a covariate, as shown in @fig-overlayed_scans below.

```{figure} #overlayed_scans
:name: fig-overlayed_scans
:width: 500px
:align: center

Overlayed Scans on Three Models
```

```{list-table}
:header-rows: 1
:name: marker-table

* - Molecular Marker 
  - Dominant/Recessive
  - Threshold 
  - Phenotype
* - D6Mit46 
  - Recessive 
  - 0.05 
  - High Gall Count
* - D8Mit88 
  - Recessive 
  - 0.05 
  - High Gall Count
* - D2Mit94 
  - Dominant 
  - 0.05 
  - High Gall Score
* - D8Mit88 
  - Recessive 
  - 0.05 
  - High Gall Weight
* - D6Mit62 
  - Dominant 
  - 0.05 
  - High Cholesterol Aggregate Crystal Formation
* - D18Mit64 
  - Recessive 
  - 0.05 
  - High Susceptibility (Augmented Data)
* - D18Mit64 
  - Recessive 
  - 0.05 
  - High Susceptibility (Original Data)
```

As we can see, the augmented data narrows down the peaks into only a few spots, showing how it counters overfitting due to the small nature of the original data set. As the original data set was much smaller than the augmented one, we can hypothesize that there are fewer peaks on the augmented model as it eliminates peaks on the original which could be due to overfitting.

## Discussion

In this study, we identified a great many locis of interest to investigate as shown in the table above.

These results are statistically significant as each of them passes the threshold value of 0.05 meaning that there is a 95 percent confidence rate for a correlation between that loci and the associated genotypes. In targeting this disease from a genomic standpoint, it may be worthy to first target those markers that are dominant for the associated phenotype. These would be D2Mit94 and D6Mit62. Research may be further made into chromosomes 2, 6, 8, and 18 as those are the most prominent chromosomes which correlate to increased gallstone susceptibility.

## Conclusion

Through our comprehensive QTL analysis exploring the genomic correlation for gallstone susceptibility, we uncovered three chromosomes of interest and 2 molecular markers of interest to target first. We conducted a robust exploration using the R/qtl package, leveraging the power of statistical methods like the Expectation-Maximization model.

Our findings illuminated several key loci with high logarithm of the odds (LOD) scores, providing significant insights into the genetic underpinnings of gallstone susceptibility. Notably, we identified loci associated with gallstone count, severity, weight, and cholesterol monohydrate crystals. The allelic variations at these loci demonstrated correlations with distinct phenotypic presentations, unraveling the complexity of genetic influences on gallstone-related traits. Furthermore, by employing an augmented dataset and comparing results with the original dataset, we sought to enhance the robustness of our analysis. The augmented data, with its increased sample size, presented a more comprehensive view of the genomic landscape associated with gallstone susceptibility. The overlay of scans from the original and augmented datasets, along with the inclusion of a covariate, provided a nuanced understanding of the genetic factors at play.

Our study contributes to the fundamental understanding of gallstone susceptibility and lays the foundation for personalized preventive strategies and therapeutic interventions. The identified QTLs harbor candidate genes that may play pivotal roles in gallstone formation, paving the way for further targeted research.

+++{"part":"data_availability"}
The data for this work was obtained from https://phenome.jax.org/projects/Lyons1.
+++

+++{"part":"acknowledgments"}
Special thanks are extended to Mr. Robert Gotwals of the North Carolina School of Science and Math and the Mouse Phenotype Database at the Jackson Laboratory.
+++

+++{"part":"author_contributions"}
This paper is solely the work of the author. All references are included in the bibliography and are cited appropriately.
+++

+++{"part":"competing-interest"}
The authors declare that they have no competing interests.
+++