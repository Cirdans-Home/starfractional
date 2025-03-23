# Routines for the solution of the Schroedinger equation

| **Routine** | **Description** |
|-------------|-----------------|
|`[Q,H] = arnoldi(A,q1,m)` | Carries out $M$ iterations of the Arnoldi iteration with $N \times N$ matrix $A$ and starting vector $\mathbf{q}$ (which need not have unit 2-norm)|
|`gen_schur_hside(m,string)` | Precomputes the representation in the Legendre basis which are used in the outer solution routines |
| `[F, k] = genCoeffMatrix_SP_interval(fV,M,I,trunc)` | Computes the coefficient matrix of $f(t)\Theta(t-s)$, with a smooth function $f$ and the Heaviside $\Theta$, in a basis of orthonormal Legendre polynomials |


> [!IMPORTANT]
> The routines `flmm2.m` and `fde12.m` are routine made by R. Garrappa
> who have been slightly adapted to work with complex dynamics by removing
> few calls to the `real()` function which in the original code was used 
> to clean spurious immaginary parts after calls to the `fft()` function.
> If you need to work with **real problems** please recover the original
> implementations. In any case, please always cite the related works:
> ```bibtex
> @article{garrappa2018numerical,
>	title={Numerical solution of fractional differential equations: A survey and a software tutorial},
>	author={Garrappa, Roberto},
>	journal={Mathematics},
>	volume={6},
>	number={2},
>	pages={16},
>	year={2018},
>	publisher={Multidisciplinary Digital Publishing Institute}
>}
>@article {MR3327641,
>    AUTHOR = {Garrappa, Roberto},
>     TITLE = {Trapezoidal methods for fractional differential equations:
>              theoretical and computational aspects},
>   JOURNAL = {Math. Comput. Simulation},
>  FJOURNAL = {Mathematics and Computers in Simulation},
>    VOLUME = {110},
>      YEAR = {2015},
>     PAGES = {96--112},
>      ISSN = {0378-4754,1872-7166},
>   MRCLASS = {65L05 (34A08 65L20)},
>  MRNUMBER = {3327641},
> MRREVIEWER = {Justin\ Steven C. Prentice},
>       DOI = {10.1016/j.matcom.2013.09.012},
>       URL = {https://doi.org/10.1016/j.matcom.2013.09.012},
>}
>```