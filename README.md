# DREAM Lambda Benchmark
This worklflow pipeline simulates biological sequencing data to test [Iota](https://github.com/swewagner/iota) as a prefilter for [Lambda](https://github.com/seqan/lambda) (version 3).

## download and installation

You can download the source code directly from GitHub:
```
git clone --recurse-submodules https://github.com/swewagner/DREAM_lambda_benchmark.git
```

Then the executables can be build:

```
mkdir build
cd build
cmake ..
make -j4
mkdir iota
```
You then need to place the `iota` binary in the `build/iota` directory and are ready to go

## usage instructions
The workflow is executed via snakemake. Each run depends on different parameters which can be set in the config yaml.

## requirements

- snakemake
- conda
- zzlib

## copyright

This workflow is free software: you can redistribute it and/or modify it under the terms of the MIT License. It is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the file [LICENSE.md](LICENSE.md) for a full text of the licence and the rights and obligations implied.

For the submodules, additional license terms apply. See the respective files in the submodules folder.

