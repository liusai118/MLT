# AI-Driven Discovery of Dual Anti-Aging and Anti-Alzheimer’s Therapeutics via PROTAC Target Deconvolution of a Super-Enhancer-Regulated Axis

> AI Discovery of Anti-Aging and Anti-AD Therapeutics

---

## Abstract

The lack of safe, durable therapeutics that act against both biological aging and Alzheimer’s disease (AD) is an unmet clinical need. To bridge this gap, we devised an AI-enabled approach that pairs rapid compound triage with mechanistic target deconvolution. Our AI-driven screening highlighted melatonin (MLT) as a promising candidate. Serum profiling of 161 human individuals confirmed an age-related fall in circulating MLT level, while subsequent in vivo and in vitro experiments showed that MLT rescues cognition, suppresses neuro-inflammation, and alleviates senescence phenotypes. PROTAC-guided chemoproteomic deconvolution next pinpointed the histone acetyltransferase p300 as MLT’s target. Integrated CUT&Tag, single-cell RNA-seq and spatial transcriptomics revealed that MLT-bound p300 cooperates with SP1 at a BMAL1 super-enhancer, elevating H3K27ac and re-engaging a circadian–epigenetic programme that links redox resilience to neuroprotection. By combining AI-driven discovery with PROTAC-based target mapping and super-enhancer-centric mechanistic resolution, our study identifies melatonin as a dual-action candidate and sets out a reproducible “AI-to-clinic” paradigm for multi-target drug innovation in aging-related neurodegeneration.

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

