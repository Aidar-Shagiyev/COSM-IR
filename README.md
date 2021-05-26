# COSM-IR
Scripts for the upcoming expansion for the [Coarse-grained Origami Structures Modeling](https://vsb.fbb.msu.ru/cosm/) web service. The expansion will feature simulating ionizing radiation-induced damage to a DNA origami structure.

## Getting Started
### Prerequisites
  * Python 3.8+
  * numpy 1.18.2+
### Using the script
For now all the arguments are hardcoded, so to transform the structure just execute cosm2cosmir.py:
```
python cosm2cosmir.py
```

## File descriptions
### Scripts
* cosm2cosmir.py: the main script; transforms COSM structure (in pdb format) into the COSM-IR structure (with simulated damage; also in pdb format).
* damage_probabilities.py: supplementary script; calculates the necessary damage probabilities.
### Examples
* example/inp: example input files for the COSM web-service
* example/out: example output files of the COSM web-service
* example/example_ir.pdb: example output file of cosm2cosmir.py
