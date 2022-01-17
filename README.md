# NMRSpecifyRegions

[![Build Status](https://github.com/RoyCCWang/NMRSpecifyRegions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/RoyCCWang/NMRSpecifyRegions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/RoyCCWang/NMRSpecifyRegions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/RoyCCWang/NMRSpecifyRegions.jl)

Given a set of frequency locations (e.g., the frequency locations from a Discrete Fourier transform of the NMR FID data) and a list of metabolites, this package filters the set such that the locations near a resonance component remains. This package relies on data from the GISSMO database for the resonance component location prediction.

See the scripts in the './examples/' folder.
