# EnGRNT

EnGRNT is a supervised ensemble learning framework for Gene Regulatory Network (GRN) inference.  
This method integrates gene expression data with topological features extracted from preliminary GRNs to train ensemble classifiers and enhance the accuracy of the inferred network structure.

This repository contains the implementation developed as part of my Master's thesis research.

---

## Published Paper

The full methodological details are available in our published paper:

Khojasteh, H., Khanteymoori, A., & Olyaee, M. H. (2021).  
**EnGRNT: Inference of gene regulatory networks using ensemble methods and topological feature extraction.**  
Informatics in Medicine Unlocked, 27, 100773.

🔗 https://www.sciencedirect.com/science/article/pii/S235291482100246X

---

## Method Overview

EnGRNT:
- Uses gene expression data
- Extracts topological features from initial GRNs
- Trains ensemble classifiers
- Improves GRN reconstruction performance

The method was evaluated on:
- Simulated datasets
- Stable-state data of *Escherichia coli (E. coli)*

With minor modifications, the framework can be adapted to other datasets.

---

## Requirements

- R version 1.1.453 or higher

---

## Academic Context

This implementation was developed during my Master's thesis and accompanies the published research work.

---

## Contact

If you have any questions or would like to collaborate, please feel free to contact:

Hakimeh Khojasteh  
khojasteh.hb@gmail.com

---

## Citation

If you use this code in your research, please cite:

```
Khojasteh, H., Khanteymoori, A., & Olyaee, M. H. (2021).
EnGRNT: Inference of gene regulatory networks using ensemble methods 
and topological feature extraction.
Informatics in Medicine Unlocked, 27, 100773.
```
