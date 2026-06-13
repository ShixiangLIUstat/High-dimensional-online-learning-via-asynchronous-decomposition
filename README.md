# High-dimensional-online-learning-via-asynchronous-decomposition
R Code of the manuscript "High-dimensional online learning via asynchronous decomposition", available at https://arxiv.org/abs/2603.20696

This repository contains R scripts in the "Code" directory that demonstrate simulation and real-data analysis.

## Files

* `R/generate.R` — Generate multiple precision matrices with similar structure for simulation studies, support Erdos-Renyi, band, block hub, and block scale-free type networks.
* `R/JGL.R` — Implementations of competing joint estimation methods: BIC-based GGL, JEM, and FJEM.
* `R/JGL MIGHT.R` — Implementation of the proposed MIGHT algorithm.
* `R/separate.R` — Implementations of competing separate estimation methods: BIC-based separate Glasso and separate nodewise regression.
* `R/simuL.R` — Code for a single simulation run (one repetition).
* `R/MainCode.R` — Main script to run MIGHT simulations: runs experiments, assembles result tables, and produces figures.
* `R/Additional simulation.R` — Additional script to run MIGHT simulations (appeared in the supplementary material): varying (n,K) in a larger p=200, with varying structure, and also heavy-tailed distributions


## Requirements
R >= 4.0
Required packages: glmnet, snowfall, Matrix, mccr, stringr, ggplot2


