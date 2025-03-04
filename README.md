# Efficient two-stage designs for benchmark dose analysis

This repository accompanies our manuscript on two stage designs and contains all code needed to replicate our results. 
Code for finding optimal designs for the second stage can be found in the design folder.
For an example of how to use the code to find designs, see the testing folder.
Simulation code can be found in the simulation folder.
Code for generating the figures can found in the figures folder.
Note that some paths are hard coded and it is recommend to change these.

This code depends on a few packages that you will likely need to install. 
Some such as ggplot2 and metaheuristicOpt are CRAN packages.
Two key packages, [ToxicR](https://github.com/NIEHS/ToxicR) and [simChef](https://yu-group.github.io/simChef/), are currently not on CRAN at the time of writing.
These packages need to be manually installed from their respective Github pages.

We also offer a Shiny app for finding the optimal second stage designs introduced in the manuscript.
A live version can be found [here]() and the source code can be found [here]().
