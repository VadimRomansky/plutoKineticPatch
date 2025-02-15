# plutoKineticPatch
module for solving kinetic convective-diffusive equation for PLUTO-4.4-2

list of changes:

1) solvig convective-diffusive equation  dF(r, p)/dt = - div D grad F(r,p) + div (u F(r,p)) + div u p/3 dF(r,p)/dp
2) template for turbulent field W(r, k)
3) Monte-Carlo particles
4) sparce output (every Nth point) (works only for dbl output)
