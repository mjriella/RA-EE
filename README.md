# RA-EE
Two Reynolds-Averaged Euler-Euler models. One is segregated and one is coupled. Additional turbulence modelling is also provided. 

Tutorials/validation cases for each chapter have been provided. Please note that at the moment the v2f model does 
not work. It was ported from a different version of OpenFOAM - if i get round to it I will fix it.

The coupled solver requires a 7x7 block therefore the VectorN for size 7 needs to be defined. Copy and paste the src inside forCompiliationCoupled into the root of your foam-extend-4.0 and recompile. There is definitely a more intelligent way of doing this but needs must.

# References
Riella, M. (2019).
”Turbulence modelling of fluid-particle interaction“
University of Exeter, PhD Thesis.

Riella, M., Kahraman, R., Tabor, G. (2018).
”Reynolds-Averaged Two-Fluid Model prediction of moderately dilute fluid-
particle flow over a backward-facing step“
Int. J. Multiphase Flow, 106:95-108. https://www.sciencedirect.com/science/article/pii/S0301932217309850

Riella, M., Kahraman, R., Tabor, G. (2019).
”Near-wall modelling in Eulerian-Eulerian Simulations“
Computers and Fluids.

Riella, M., Kahraman, R., Tabor, G. (2019).
”Inhomogeneity and anisotropy in Eulerian-Eulerian near-wall modelling“
Int. J. Multiphase Flow, 114:9-18. https://www.sciencedirect.com/science/article/pii/S0301932218304932

Riella, M., Kahraman, R., Tabor, G. (2019).
”Fully-coupled pressure-based two-fluid solver for the solution of turbulent
fluid-particle systems“
Computer and Fluids.
