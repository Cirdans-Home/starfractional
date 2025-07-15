# ⋆-Product Approach for Analytical and Discrete Solutions of Nonautonomous Linear Fractional Differential Equations

This repository contains the accompanying code for the article: "⋆-Product Approach for Analytical and Discrete Solutions of Nonautonomous Linear Fractional Differential Equations"

## Collaborators

- Fabio Durastante [:email:](mailto:fabio.durastante@unipi.it)
- Pierre-Louis Giscard [:email:](mailto:giscard@univ-littoral.fr)
- Stefano Pozza [:email:](mailto:pozza@karlin.mff.cuni.cz)

## External codes

To run the code in this repository you need the `Chebfun` code suite. This can be obtained by running the following commands in the MATLAB command window:
```matlab
unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath
```

> [!IMPORTANT]
> If you use this code in a scientific publication it is important to cite:
> ```
> T. A. Driscoll, N. Hale, and L. N. Trefethen, editors, Chebfun Guide, Pafnuty Publications, Oxford, 2014. 
> ```

The other piece of external code which is used are the Fractional Linear
Multistep Method implemented by R. Garrappa, which are used as comparison
with the procedure presented here.

> [!IMPORTANT]
> If you use such codes in a scientific publication it is important to cite:
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
