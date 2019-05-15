# SuCOS

SuCOS: a RDKit-based shape overlap metric. Combines volume overlap and featuremap scores to give a combined overlap score similar to OpenEye's Tversky Combo.

For calculation of chemical feature overlap, most of the code is based on this RDKit blog post about [using feature maps](http://rdkit.blogspot.com/2017/11/using-feature-maps.html)

If you use this code, please cite our ChemRxiv [manuscript](https://chemrxiv.org/articles/SuCOS_is_Better_than_RMSD_for_Evaluating_Fragment_Elaboration_and_Docking_Poses/8100203).

## Getting Started

### Prerequisites

* [RDKit](http://www.rdkit.org/) 
* Numpy 

## Running the tests

Run the unit tests.

```
> python test_SuCOS.py
```

## Command line help.
```
> python calc_SuCOS.py -h

usage: calc_SuCOS.py [-h] [--lig1 LIG1] [--lig2 LIG2] [--write] [--return_all]
                     [--score_mode {all,closest,best}]

run SuCOS

optional arguments:
  -h, --help            show this help message and exit
  --lig1 LIG1           first ligand, in sdf or mol2 file format.
  --lig2 LIG2           second ligand(s), in sdf or .sdf.gz file format.
  --write               writes the SuCOS score into a sdf file with suffix
                        _SuCOS_score.sdf
  --return_all
  --score_mode {all,closest,best}
                        choose the scoring mode for the feature map, default
                        is all.

```
## An example
```
> python calc_SuCOS.py --lig1 test_data/4e3g_lig.sdf --lig2 test_data/benzene.sdf 

********************************
SuCOS score:	0.843867
********************************
```
