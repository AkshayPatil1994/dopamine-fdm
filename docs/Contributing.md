# Contributing & Citing

← [[Home|Home]]

## License

This program is free software: you can redistribute it and/or modify it under the terms
of the GNU Affero General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

See [LICENSE](../LICENSE) for the full text.

## Reporting issues / contributing

- If you spot a mistake anywhere in the docs (including this wiki), please open an issue
  on GitHub.
- Bug reports and pull requests are welcome; there is no separate CONTRIBUTING template
  at present, so a descriptive issue or PR is the way to start.

## Acknowledgments

The original form of the channel flow solver was shared by Adrian Lozano-Duran during
the author's PhD; his generous help and support with the code are gratefully
acknowledged.

## Citing this solver

If you use `dopamine-fdm` in your research, please consider citing the following
publications (also machine-readable in [`CITATION.cff`](../CITATION.cff)):

1. Lozano-Durán, A., & Bae, H. J. (2019). Characteristic scales of Townsend's
   wall-attached eddies. *Journal of Fluid Mechanics*, 868, 698–725.
   doi:[10.1017/jfm.2019.209](https://doi.org/10.1017/jfm.2019.209) — *original DNS
   solver*
2. Patil, A., & Fringer, O. (2022). Drag enhancement by the addition of weak waves to a
   wave-current boundary layer over bumpy walls. *Journal of Fluid Mechanics*, 947, A3.
   doi:[10.1017/jfm.2022.628](https://doi.org/10.1017/jfm.2022.628) — *IBM + DNS*
3. Patil, A., & García-Sánchez, C. (2025). Should we care about the spatial
   heterogeneity in coral reefs under unidirectional turbulent flows?
   arXiv:[2506.03021](https://arxiv.org/abs/2506.03021) [physics.flu-dyn] — *improved
   IBM*
4. Patil, A., Paranjothi, U. C. K., & García-Sánchez, C. (2025). GenSDF: An MPI-Fortran
   based signed-distance-field generator for computational fluid dynamics applications.
   *SoftwareX*, 30, 102117. — *GenSDF* (see
   [[Pre- and Post-Processing Tools § GenSDF|Tools#gensdf]])

<details>
<summary>BibTeX</summary>

```bibtex
@article{lozanoduran2019characteristic,
  title={Characteristic scales of {T}ownsend's wall-attached eddies},
  author={Lozano-Dur{\'a}n, Adri{\'a}n and Bae, Hyunji Jane},
  journal={Journal of Fluid Mechanics},
  volume={868},
  pages={698--725},
  year={2019},
  publisher={Cambridge University Press},
  doi={10.1017/jfm.2019.209}
}

@article{patil2022drag,
  title={Drag enhancement by the addition of weak waves to a wave-current boundary layer over bumpy walls},
  author={Patil, Akshay and Fringer, Oliver},
  journal={Journal of Fluid Mechanics},
  volume={947},
  pages={A3},
  year={2022},
  publisher={Cambridge University Press},
  doi={10.1017/jfm.2022.628}
}

@misc{patil2025carespatialheterogeneitycoral,
  title={Should we care about the spatial heterogeneity in coral reefs under unidirectional turbulent flows?},
  author={Patil, Akshay and Garc{\'i}a-S{\'a}nchez, Clara},
  year={2025},
  eprint={2506.03021},
  archivePrefix={arXiv},
  primaryClass={physics.flu-dyn},
  url={https://arxiv.org/abs/2506.03021}
}

@article{patil2025gensdf,
  title={{GenSDF}: An {MPI-Fortran} based signed-distance-field generator for computational fluid dynamics applications},
  author={Patil, Akshay and Paranjothi, Uma Chandrika Karrothu and Garc{\'i}a-S{\'a}nchez, Clara},
  journal={SoftwareX},
  volume={30},
  pages={102117},
  year={2025},
  publisher={Elsevier}
}
```

</details>

## Methods reference table

The methods implemented in this solver draw on the following published works. Please
cite the relevant papers when publishing results obtained with `fdm-dopamine`.

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

See the fuller in-context references throughout [[Numerics & Governing Equations|Numerics]].
