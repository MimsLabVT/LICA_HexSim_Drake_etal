Combining graph theory and spatially-explicit, individual-based models to improve
invasive species control strategies at a regional scale
================

## Table of Contents

-   Introduction
-   Code for generating simulated landscapes
-   Code for combining census files
-   Code for analysis of derived HexSim output files

## Introduction and Caveats

This markdown includes code used to collate, shape, and analyze HexSim
output data associated with the manuscript: Combining graph theory and spatially-explicit,
individual-based models to improve invasive species control strategies at a regional scale
This manuscript focuses on exploring optimal management strategies of the American bullfrog, *Lithobates catesbaeinus*,
informed by spatial data from the Huachuca Mountains and Canelo Hills in
southeastern Arizona.

The first section details how we simulated multiple landscapes for use in the IBM.
As a result of the need of multiple replicates across landscape
iterations and for each sensitivity analysis and for each simulated
management scenario, there was a bounty of data needed to be shaped and
combined to allow for downstream analyses in the final sections. In the second section, we
present the code to replicate the analyses reported in the manuscript.


While we present these code for transparency of methodology and
potential replication, newer version of HexSim have added increased
functionality rendering some of the workarounds (which were required the
carpentry of the data we analyzed derived from the individually based
model and its iteration) not necessarily needed. Read the annotations to
better understand different aspects of the code and what was performed.
We would recommend these scripts being embedded into an R project for
convenience of directories working to your benefit.

Preferred Citation:

Drake, J.C., O'Malley, G., Kraft, J., Mims, M.C. 2024. Combining graph theory and spatially-explicit,
individual-based models to improve invasive species control strategies at a regional scale. *Landscape Ecology* (In Revision).

##

