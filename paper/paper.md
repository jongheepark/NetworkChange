---
title: 'NetworkChange: Analyzing Network Changes in R'
tags:
  - R
  - longitudinal networks
  - latent space 
  - hidden Markov models
  - Bayesian inference
authors:
- affiliation: 1
  name: Jong Hee Park
  orcid: 0000-0002-4598-3456
- affiliation: 2
  name: Yunkyu Sohn
  orcid: 0000-0002-9496-1907
affiliations:
- index: 1
  name: Department of Political Science and International Relations, Seoul National University
- index: 2
  name: Faculty of Political Science and Economics, Waseda University
date: 9 July 2020
bibliography: paper.bib
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Journal of Open Source Software.
---

# Summary

**NetworkChange** is an R package that detects structural
changes in longitudinal network data using the latent space approach [@hoff2002latent].
Based on the Bayesian multi-array representation of longitudinal
networks [@Hoff2015], **NetworkChange** performs
Bayesian hidden Markov analysis to discover changes in structural
network features across temporal layers using the hidden Markov model formulation. **NetworkChange** can detect
various forms of changes in network structure such as block-splitting, block-merging,
and core-periphery changes. **NetworkChange** also provides functions
for model diagnostics using WAIC, average loss, and log marginal
likelihoods as well as visualization tools for dynamic analysis results
of longitudinal networks. 

# Statement of need 

The package is designed for R users who need to analyize longitudinal network data to discover latent node-level characteristics including cases when there are discrete changes of the underlying states governing the node-level characteristics. In addition to functions for the statistical analysis, **NetworkChange** provides visualization functions for summary of the analysis results. The complete guide for using core functions of the package is presented at https://github.com/jongheepark/NetworkChange as its vignette with an empirical data set analysis example. Methodological details of the package were published as an article at *Bayesian Analysis* [@ParkSohn2020].

# Acknowledgements

This work was supported by the Japan Society for the Promotion of Science Early-Career Scientists Grant [JP19K13606] to Y.S.

# References