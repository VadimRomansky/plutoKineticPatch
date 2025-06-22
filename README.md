# plutoKineticPatch
module for solving kinetic convective-diffusive equation for PLUTO-4.4

list of changes:

1) solvig convective-diffusive equation  dF(r, p)/dt = - div D grad F(r,p) + div (u F(r,p)) + div u p/3 dF(r,p)/dp
2) template for turbulent field W(r, k)
3) Monte-Carlo particles
4) sparce output (every Nth point) (works only for dbl output)
5) renaming SZ and SZ
6) fixed bug with reading multiple files in pypluto and added some plots

Install

merge everything in PLUTO directory
