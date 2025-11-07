<!-- Title -->
<h1 align="center">
  FLOW2
</h1>

<!-- Information badges -->
<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg" />
  </a>
  <a href="https://fortran-lang.org">
    <img alt="Fortran" src="https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat" />
  </a>
</p>

<p align="center">
  <img width="300" src="img/nstx_flowcurves.png" alt="NSTX flow curves">
</p>

FLOW2 was developed at the University of Rochester. FLOW2 is written primarily in Fortran and is used to study the equilibrium properties of toroidal devices, such as tokamaks, in conditions relevant to present day experiments. 

The most unique feature of the code is the ability to study flow-dependent equilibria with distinction between ion and electron properties. The two-fluid model differentiates the fluid properties of ions and electrons, in particular pressures, velocities, and in principle densities. Since electron inertia is negligible with respect to ion inertia, plasma macroscopic mass flow is determined by ion velocities. An effect of this fact is that, contrary to the MHD case, plasma flow will be on surfaces that are **not** magnetic flux surfaces. As a consequence, there will in general be a finite plasma normal velocity to any magnetic surface, a fact which is not allowed in MHD theory. [R. Betti, J. P. Freidberg, PoP 7, 2439 (2000)](https://pubs.aip.org/aip/pop/article-abstract/7/6/2439/103410/Radial-discontinuities-in-tokamak?redirectedFrom=fulltext).

Despite the obvious relevance of the study of equilibrium in the presence of strong flow, little work has been done on the subject. The main macroscopic difference between MHD and two-fluid equilibria is the presence of a finite plasma velocity *normal* to magnetic surfaces. The analysis of the effect of this new (with respect to MHD) element introduced by the two-fluid model is the study of ongoing research.

The code input is formulated in terms of “quasi-physical” free functions, which allow for an intuitive and user-friendly interface. For a more detailed description of the code and its inputs, see [L. Guazzotto, R. Betti, PoP 22, 092503 (2015)](https://pubs.aip.org/aip/pop/article/22/9/092503/109618/Two-fluid-equilibrium-with-flow-FLOW2).

## Links

* Project homepage: <https://www.auburn.edu/cosam/faculty/physics/guazzotto/research/FLOW2_main.html>
* Author homepage: <https://www.auburn.edu/cosam/faculty/physics/guazzotto/research/>
* Code repository: <https://github.com/l-guazzotto/FLOW2>