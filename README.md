# Precision-Recall Assessment of Protein Function Prediction

## Introduction
Critical Assessment of Function Annotation (CAFA), is a community-wide challenge designed to provide a large-scale assessment of computational methods dedicated to predicting protein function.

More information can be found at http://biofunctionprediction.org/cafa/ as well as the CAFA2 paper (Jiang et al, 2016)

This toolset provides an assessment for CAFA submissions based on precision and recall. 


## Dependencies
 - Python 2.7 or Python 3
 - Python packages:
 1. Biopython (Bio)
 2. yaml
 3. matplotlib 
 4. seaborn

## Main Functions
 - `assess_main.py` 
 - `plot.py`


## Auxiliary Functions 
 CAFA3 released its [protein targets](https://www.synapse.org/#!Synapse:syn6172284) in September 2016. Each protein target has a unique CAFA3 ID. To run the above assessment function, each protein should be represented by its CAFA3 ID. However, the benchmark proteins generated by the [benchmark creation tool](https://github.com/nguyenngochuy91/CAFA_benchmark) are identified by UniProt Accession IDs. 
Therefore, we here provide functions to convert between UniProt IDs and CAFA3 IDs. We also provide a function that converts benchmark files generated by the [benchmark creation tool](https://github.com/nguyenngochuy91/CAFA_benchmark) to a benchmark folder that can feed into this program.
 - `benchmark_folder.py`

## Examples


## References
Jiang, Yuxiang, et al. "An expanded evaluation of protein function prediction methods shows an improvement in accuracy." Genome biology 17.1 (2016): 184.

http://biofunctionprediction.org/cafa/
