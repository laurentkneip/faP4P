# Project Name

faP4P

## Table of Contents

- [Brief](#brief)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Brief

This is an implementation of the uncalibrated P4P algorithm of Yang Guo published in 2013 in JMIV. The algorithm solves for pose, unknown focal length, and aspect ratio. 

Paper:

A Novel Solution to the P4P Problem for an Uncalibrated Camera, Yang Guo, J Math Imaging Vis (2013) 45:186–198, DOI 10.1007/s10851-012-0360-0


## Installation

```bash
# Example installation steps or commands
git clone https://github.com/laurentkneip/fap4p.git
cd fap4p
mkdir build
cd build
cmake ..
make
```

Note that the code depends on Eigen!

## Usage

Just go to the bin folder and type ./test . A random test is generated and the solutions are printed into the console every time the program is run.

## License

GNU GENERAL PUBLIC LICENSE