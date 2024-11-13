# VPalm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PalmStudio.github.io/VPalm.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PalmStudio.github.io/VPalm.jl/dev/)
[![Build Status](https://github.com/PalmStudio/VPalm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PalmStudio/VPalm.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/PalmStudio/VPalm.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PalmStudio/VPalm.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

VPalm is an automaton that builds 3d mockups of palm plants from architectural parameters and allometric equations.

The model also integrates a biomechanical model to compute the leaf bending and torsion using the biomass of each leaf.


To do:

- [ ] Change the way we compute the height of the internodes, atm it needs to know the total number of internodes to compute the height of each internode, but this changes every time there is a new internode, which is false because we don't have secondary growth for palm trees.