# SuCOS

SuCOS is an RDKit-based overlap measure that combines volumetric shape and pharmacophoric features to give a combined overlap score similar to OpenEye's Tversky Combo and [Malhotra and Karanicolas' COS measure](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.6b00725). 

For calculation of chemical feature overlap, most of the code is based on this RDKit blog post about [using feature maps](http://rdkit.blogspot.com/2017/11/using-feature-maps.html)

If you use this code, please cite our ChemRxiv [manuscript](https://chemrxiv.org/articles/SuCOS_is_Better_than_RMSD_for_Evaluating_Fragment_Elaboration_and_Docking_Poses/8100203), where we show that SuCOS is better than RMSD and PLIF similarity for evaluating binding mode conservation in fragment-based drug discovery. 

## Getting Started

### Prerequisites

* [RDKit](http://www.rdkit.org/) 
* Numpy 

## Running the tests

You can run the unit tests by typing:

```
> python test_SuCOS.py
```

## Command line help.

To list the options available in calc_SuCOS.py, type:

```
> python calc_SuCOS.py -h

usage: calc_SuCOS.py [-h] [--lig1 LIG1] [--lig2 LIG2] [--write] [--return_all]
                     [--score_mode {all,closest,best}]

run SuCOS

optional arguments:
  -h, --help            show this help message and exit
  --lig1 LIG1           the smaller/reference ligand, in sdf or mol2 file
                        format.
  --lig2 LIG2           the larger/query ligand(s), in sdf or .sdf.gz file
                        format.
  --write               writes the SuCOS score into a sdf file with suffix
                        _SuCOS_score.sdf
  --return_all
  --score_mode {all,closest,best}
                        choose the scoring mode for the feature map, default
                        is all.

```
## Example

To find the overlap between 4e3g_lig.sdf and benzene.sdf in the test_data directory:

```
> python calc_SuCOS.py --lig1 test_data/4e3g_lig.sdf --lig2 test_data/benzene.sdf 

********************************
SuCOS score:	0.843867
********************************
```
## License

SuCOS is licensed under the MIT license.
