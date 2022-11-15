# *-fractional: solving linear nonautonomous fractional differential equations
The *-Lanczos method for Fractional Ordinary Differential Equations

## Getting the code

To obtain the code run
```
git clone --recurse-submodules git@github.com:Cirdans-Home/starfractional.git
```
- The `--recurse-submodules` options downloads also the GIST in the `pynotebook` folder. In any case, the folder also contains the clean sources of the notebooks. So unless you need to view them online, you can make a simple clone and open them in your local environment. In this second case, the environment must be configured to use [SageMath](https://www.sagemath.org/).
- To run the `Matlab` code you need some of functionalities available in [Chebfun](https://www.chebfun.org/), so you need to have it installed. It can be done easily by running in the `Matlab` command-line:
```
unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath
```
