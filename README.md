# COSM-IR
Scripts for the upcoming expansion for the [Coarse-grained Origami Structures Modeling](https://vsb.fbb.msu.ru/cosm/) web service. The expansion will feature simulating ionizing radiation-induced damage to a DNA origami structure.

## Getting Started
### Prerequisites
  * Python 3.8+
  * numpy 1.18.2+
### Using the script
To transform the structure execute cosm2cosmir.py:
```
python cosm2cosmir.py -i example/out/example.pdb -o example/example_ir.pdb -r example/out/example.r -ro example/example_ir.r -s example/inp/example_scaffold.txt -m example/out/example.map -d 100
```

## File descriptions
### Force Field
* cosmir.ff contains modified COSM force field with new particles for modeling damaged base pairs.
### Scripts
* cosm2cosmir.py: the main script; transforms COSM structure (in pdb format) into the COSM-IR structure (with simulated damage; also in pdb format).
* damage_probabilities.py: supplementary script; calculates the necessary damage probabilities.
### Examples
* example/inp: example input files for the COSM web-service
* example/out: example output files of the COSM web-service
* example/example_ir.pdb: example output file of cosm2cosmir.py
