# WeaverLangmuir

WL1989.for	First Commit	6 months ago
wl1989kdcalc.py	First Commit	6 months ago
wl1989models.py	First Commit	6 months ago
wl1989stoich.py	First Commit	6 months ago
wlState.py

This folder is an updated and modified version of Weaver and Langmuir 1989 published in Computers & Geosciences. The original code is WL1989.for.

I rewrote the code in python and updated the Kd calculations. The set up of the code is generally organized the same as the original Fortran Code.

wl1989models.py is the code that calculates crystallization of a magma chamber. This code calls wlState which calls wl1989kdcalc and wl1989stoich.
