# References

The methods implemented in this solver draw on the following published works.  Please cite the relevant papers when publishing results obtained with `fdm-dopamine`.

| Method | Citation | DOI |
|--------|----------|-----|
| Time integration (RK3) | Wray, A.A. (1990). *Minimal storage time advancement schemes for spectral methods*. NASA Ames Report. | — |
| Fractional-step projection | Kim, J. & Moin, P. (1985). J. Comput. Phys. **59**, 308–323. | [10.1016/0021-9991(85)90148-2](https://doi.org/10.1016/0021-9991(85)90148-2) |
| Spectral Poisson solver (FFTW3) | Frigo, M. & Johnson, S.G. (2005). Proc. IEEE **93**(2), 216–231. | [10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301) |
| MPI pencil decomposition ([2decomp&fft](https://github.com/2decomp-fft/2decomp-fft)) | Li, N. & Laizet, S. (2010). *2DECOMP&FFT – A Highly Scalable 2D Decomposition Library for FFT-based Simulations*. Cray User Group 2010. | — |
| Ghost-cell IBM | Tseng, Y.-H. & Ferziger, J.H. (2003). J. Comput. Phys. **192**(2), 593–623. | [doi:10.1016/j.jcp.2003.07.024](https://doi.org/10.1016/j.jcp.2003.07.024) |
| Vreman SGS model | Vreman, A.W. (2004). Phys. Fluids **16**(10), 3670–3681. | [10.1063/1.1785131](https://doi.org/10.1063/1.1785131) |
| van Leer MUSCL limiter (scalar) | van Leer, B. (1974). J. Comput. Phys. **14**(4), 361–370. | [10.1016/0021-9991(74)90019-9](https://doi.org/10.1016/0021-9991(74)90019-9) |
| Settling velocity (sediment) | Soulsby, R.L. (1997). *Dynamics of Marine Sands*. Thomas Telford. | ISBN 978-0-7277-2584-5 |
| Reynolds stress budget | Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press. | [10.1017/CBO9780511840531](https://doi.org/10.1017/CBO9780511840531) |
| Reichardt IC profile | Reichardt, H. (1951). Z. Angew. Math. Mech. **31**(7), 208–219. | [10.1002/zamm.19510310704](https://doi.org/10.1002/zamm.19510310704) |
| Staggered MAC grid | Harlow, F.H. & Welch, J.E. (1965). Phys. Fluids **8**(12), 2182–2189. | [10.1063/1.1761178](https://doi.org/10.1063/1.1761178) |
| PLIC-VOF, split advection | Youngs, D.L. (1982). *Time-dependent multi-material flow with large fluid distortion*. In Numerical Methods for Fluid Dynamics, Academic Press. Rider, W.J. & Kothe, D.B. (1998). J. Comput. Phys. **141**, 112–152. | [10.1006/jcph.1998.5906](https://doi.org/10.1006/jcph.1998.5906) |
| Conservative direction-split VOF (compression term) | Weymouth, G.D. & Yue, D.K.P. (2010). J. Comput. Phys. **229**(8), 2853–2865. | [10.1016/j.jcp.2009.12.018](https://doi.org/10.1016/j.jcp.2009.12.018) |
| Consistent mass–momentum transport | Rudman, M. (1998). Int. J. Numer. Meth. Fluids **28**, 357–378. Vaudor, G. et al. (2017). Comput. Fluids **152**, 204–216. Pal, S., Fuster, D. & Zaleski, S. (2021). arXiv:2101.04142. Zeng, Y. et al. (2023). J. Comput. Phys. **478** (consistent adaptive level-set framework). | [10.1016/j.compfluid.2017.04.017](https://doi.org/10.1016/j.compfluid.2017.04.017) |
| WENO5-Z | Jiang, G.-S. & Shu, C.-W. (1996). J. Comput. Phys. **126**, 202–228. Borges, R. et al. (2008). J. Comput. Phys. **227**, 3191–3211. | [10.1006/jcph.1996.0130](https://doi.org/10.1006/jcph.1996.0130), [10.1016/j.jcp.2007.11.038](https://doi.org/10.1016/j.jcp.2007.11.038) |
| Height-function curvature | Popinet, S. (2009). J. Comput. Phys. **228**, 5838–5866. | [10.1016/j.jcp.2009.04.042](https://doi.org/10.1016/j.jcp.2009.04.042) |
| Balanced-force surface tension, CSF | Brackbill, J.U., Kothe, D.B. & Zemach, C. (1992). J. Comput. Phys. **100**, 335–354. Francois, M.M. et al. (2006). J. Comput. Phys. **213**, 141–173. | [10.1016/0021-9991(92)90240-Y](https://doi.org/10.1016/0021-9991(92)90240-Y), [10.1016/j.jcp.2005.08.008](https://doi.org/10.1016/j.jcp.2005.08.008) |
| Harmonic mean of the viscosity at interface edges | Tryggvason, G., Scardovelli, R. & Zaleski, S. (2011). *Direct Numerical Simulations of Gas–Liquid Multiphase Flows*. Cambridge University Press. | [10.1017/CBO9780511975264](https://doi.org/10.1017/CBO9780511975264) |
| Reversed-vortex / deformation tests | LeVeque, R.J. (1996). SIAM J. Numer. Anal. **33**, 627–665. Enright, D. et al. (2002). J. Comput. Phys. **183**, 83–116. | [10.1137/0733033](https://doi.org/10.1137/0733033), [10.1006/jcph.2002.7166](https://doi.org/10.1006/jcph.2002.7166) |
| Wave generation, relaxation zones | Fenton, J.D. (1988). J. Waterway, Port, Coastal, Ocean Eng. **114**, 108–116 (stream function). Jacobsen, N.G., Fuhrman, D.R. & Fredsøe, J. (2012). Int. J. Numer. Meth. Fluids **70**, 1073–1088 (relaxation zones). | [10.1002/fld.2726](https://doi.org/10.1002/fld.2726) |
| Capillary-wave reference | Prosperetti, A. (1981). Phys. Fluids **24**, 1217–1223. | [10.1063/1.863522](https://doi.org/10.1063/1.863522) |
| Dam-break reference | Martin, J.C. & Moyce, W.J. (1952). Phil. Trans. R. Soc. A **244**, 312–324. Koshizuka, S. & Oka, Y. (1996). Nucl. Sci. Eng. **123**, 421–434. | [10.1098/rsta.1952.0006](https://doi.org/10.1098/rsta.1952.0006) |
| Taylor–Green vortex reference | van Rees, W.M. et al. (2011). J. Comput. Phys. **230**, 2794–2805. | [10.1016/j.jcp.2010.11.031](https://doi.org/10.1016/j.jcp.2010.11.031) |
