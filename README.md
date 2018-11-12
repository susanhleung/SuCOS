# SuCOS

SuCOS: a RDKit-based shape overlap metric. Combines volume overlap and featuremap scores to give a combined overlap score similar to OpenEye's Tanicombo.

## Getting Started

### Prerequisites

The latest RDKit. 

### Installing

TODO

## Running the tests

Run the unit tests.

```
> python test_SuCOS.py
```

## Command line help.
```
> python calc_SuCOS.py -h

usage: calc_SuCOS.py [-h] [--lig1 LIG1] [--lig2 LIG2]

run SuCOS

optional arguments:
  -h, --help   show this help message and exit
  --lig1 LIG1  first ligand
  --lig2 LIG2  second ligand
```
## An example
```
> python calc_SuCOS.py --lig1 test_data/4e3g_lig.sdf --lig2 test_data/benzene.sdf 

********************************
SuCOS score:	0.843867
********************************
```
## Built With

* [RDKit](http://www.rdkit.org/) 

## Authors

* **Susan Leung**