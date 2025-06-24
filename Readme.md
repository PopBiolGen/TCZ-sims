# Spread model for toads

This project started with the code used by Tingley *et al.* (2013, *J. Appl. Ecol* 50: 129-137).  This legacy code was not under version control.

We build off this code base to run scenario testing for the proposed Toad Containment Zone and to estimate the timing of toad arrival in the Pilbara (under a do nothing scenario), and at the head of the Toad Containment Zone.

The bulk of the code is for running the revised version of the Tingley *et al.* model as described in Dunlop *et al.* (submitted 2025).

With regard to Dunlop *et al*, that paper uses three methods for estimating the arrival time of toads in the Pilbara. Method 1 is a simple division operation and is reported in the paper.  Method 2 takes account of local rainfall and is slightly more complex.  The script for that is `src/pilbara-impact-sims/spread-rate-rainfall.R`.

Method 3 is the update of the Tingley model.  The workhorse here is `src/pprocess_functions.R`, with scripts in subfolders calling these functions against different scenarios.  The `pilbara-impact-sims` subfolder contains the simulations used to estimate the arrival time to the pilbara.
