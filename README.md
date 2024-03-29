# Serological-analysis-in-humans-in-Malaysian-Borneo-suggests-prior-exposure-to-H5-avian-influenza
Cases of highly pathogenic H5 avian influenzas are on the rise globally. The reported case fatality rate in humans is nearly 60%. Malaysian Borneo is a stopover site for wild birds migrating between Russia and Australia. Although no human cases of H5 avian influenza infection have been reported in Malaysian Borneo to date, highly pathogenic H5N1 has circulated in poultry and migratory avian species transiting through this region. Malaysian Borneo has seen recent deforestation and destruction of natural shorebird habitats, which may increase the proximity between humans and avian species. We hypothesise that higher rates of human-animal contact, caused by this habitat destruction, will increase the likelihood of potential zoonotic spillover events. In 2015, an environmentally stratified cross-sectional survey was conducted collecting geolocated questionnaire data on potential risk factors for emerging zoonotic disease in 10,100 individuals. We performed a serological survey of these individuals via ELISAs, pseudotyped neutralisation assays and a cross-reactivity depletion assay. We found evidence of H5 neutralisation that persisted following depletion of seasonal H1/H3 binding antibodies from the plasma. The presence of these antibodies suggests that some individuals living near these migratory sites may have been exposed to H5 influenza. There is a spatial overlap between individuals displaying high H5 binding and the distributions of these migratory species along with shared environmental risk factors. We have developed a novel surveillance approach including both spatial and serological data to detect potential spillover events, highlighting the urgent need to study cross-species pathogen transmission in these migratory zones.

The attached sample code file includes:

1. Spatial autocorrelation determined by Moran's I analyses
2. Multivariate analysis to identify potential environmental covariates
3. Geostatistical modelling of areas with high probability of contact with specific outcomes and shared/separate spatial effects

Data can only be shared with approval from relevant ethics committees. For further information, please contact authors.

**Instructions for use**
**System requirements:**
R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)

RStudio 2022.02.3+492 "Prairie Trillium" Release (1db809b8323ba0a87c148d16eb84efe39a8e7785, 2022-05-16) for macOS
Mozilla/5.0 (Macintosh; Intel Mac OS X 13_5_1) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36

**Installation guide**
Attached files can be downloaded and directly imported into RStudio with read.csv().

**Simulated data/demo**
An input demo dataset is attached included data on wild bird contact (real) and H5 binding (simulated) with environmental variables and simulated GPS points. A small prediction dataset is also included with simulated GPS points. Expected run time varies based on the number of covariates and whether a spatial fixed effect is included.
