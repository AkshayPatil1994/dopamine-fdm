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
