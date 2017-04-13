# Precision-Recall Assessment of Protein Function Prediction

## Introduction
Critical Assessment of Function Annotation (CAFA), is a community-wide challenge designed to provide a large-scale assessment of computational methods dedicated to predicting protein function.

The third CAFA is currently closed for submission. More information can be found at http://biofunctionprediction.org/cafa/.

This toolset provides an assessment for CAFA submissions based on precision and recall. 


## Dependencies
 - Python 2.7 and above
 - Python packages:
 1. Biopython (Bio)
 2. matplotlib
 3. seaborn
 
## Parameters
- Positional Parameter: filepath to the prediction file that you wish to assess.
Prediction file should follow CAFA format, in both filename and file content. Detailed format requirements are listed here: https://www.synapse.org/#!Synapse:syn5840147/wiki/402192.
You can submit as many files as you like in one run, they will be plotted on the same figure, as well as individually. 

- Other parameters:
1. Benchmark Type: We define two types of benchmarks in CAFA. 
No knowledge (NK) benchmarks are those proteins that have no experimental annotation in all three ontologies (BPO, CCO and MFO) at submission deadline, and gained experimental annotation in the ontology of interest.
Limited knowledge (LK) benchmarks are those proteins that have experimental annotation in one or two ontologies, but not in the one of interest at submission deadline, and gained experimental annotation in the ontology of interest.
2. Evaluation Mode: Full evaluation mode considers the entire set of benchmark proteins, while partial mode only considers a subset of the benchmark protein that has been predited by the CAFA team.
3. Benchmark Folder: Folder containing benchmark proteins and their gained experimental annotations. Default CAFA 2 benchmark folder is provided. Customized benchmark folder should follow CAFA 2 structure.
4. GO.obo File Path: Gene Ontology file used. Default is the one used for CAFA 2 evaluation.
5. Smooth: Option to have the PR curves smoothed. Recommended if plotting multiple curves on one figure.

## Execution

0. Install python
1. Download the package and cd to the main directory of the package in command-line console
2. Type `python precrec_main.py -h` for usage on this tool

## Example

`python precrec_main.py ./Doegroup_1_10116.txt -t 'type1' -m 'full' -s 'N' `
`python precrec_main.py ./Doegroup_2_9606.txt .Doegroup_1_9606.txt -t 'type1' -m full -title 'Doegroup_1_vs_2' -s Y`

## References
Jiang, Yuxiang, et al. "An expanded evaluation of protein function prediction methods shows an improvement in accuracy." Genome biology 17.1 (2016): 184.

http://biofunctionprediction.org/cafa/