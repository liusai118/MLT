# PTD-DEP: Pathway and Transcriptome-Driven Drug Efficacy Predictor

> Dual-modality AI for prioritizing neuroprotective compounds in Alzheimer’s disease (AD) and aging.

---

## Overview

We develop PTD-DEP, a dual-modality AI that couples ML-based pathway screening with deep transcriptomic modeling to prioritize neuroprotective compounds for AD and aging.  
The pathway layer (six algorithms) attains mean ROC-AUC > 0.92, while the transcriptomic layer—aligning L1000 profiles to brain-like cell types and combining a 64-D autoencoder with a two-layer GCN on SMILES—predicts expression profiles with Pearson = 0.968, enabling quantitative multi-target ranking.

---

## Results

### Result 1. Construction of deep learning-based PTD-DEP model
We develop PTD-DEP, a dual-modality AI that couples ML-based pathway screening with deep transcriptomic modeling to prioritize neuroprotective compounds for AD and aging. The pathway layer (six algorithms) attains mean ROC-AUC > 0.92, while the transcriptomic layer—aligning L1000 to brain-like cell types and combining a 64-D autoencoder with a two-layer GCN on SMILES—predicts expression profiles with Pearson = 0.968, enabling quantitative multi-target ranking. *(Fig. 1a–e)*

#### P2T-DEP — Used for training PTD-DEP  
Pipeline/scripts to train both ML and DL components.

#### ML — The ML layer of PTD-DEP  
Multi-task, multi-algorithm pathway prediction.

#### DL — The DL layer of PTD-DEP  
Autoencoder + GCN (SMILES-to-graph) transcriptomic predictor.

---

### Result 5. Melatonin improved cellular senescence in mesenchymal stem cells
Melatonin alleviated senescence phenotypes in MSCs, supporting PTD-DEP-guided prioritization. 

---

### Result 6. Super-enhancer-driven p300/SP1 potentiated **BMAL1** transcription
A super-enhancer program engaging p300/SP1 augments **BMAL1** transcription. 

- **Cut&Tag** — Used for Cut&Tag analysis  
- **MAPE_PPI** — Used for PPI prediction  
- **RNAseq** — Used for RNA-seq analysis

---

## Repository Structure

