## LOG

This version is improved from 180820 version. It has a better performance now (a little readability sacrified);
using analytical solutions of the household lifetime decision-making problems, allowing negative assets;

an extra testing script *test_Analytical.jl* is provided. It shows Euler equation always works on the same dataset and solutions starting from different ages are consistent.

## DESCRIPTION
1. updated *NumAlgo*, added *Bisection()* to do bisection searching (p.s.: in another developing version, it has been replaced by equivalent methods in *Roots* package, this package is so common that it does not bring extra dependencies)
2. added *Policy_Analytical* & *Policy_Analytical_Retired* modules; the *PolicyFunctions* module (dynamic programming based) has been deprecated (but not deleted)
3. added *New_DepartableUtilityFunction* and other attached documents which discussed mathematics, algorithms and notations in the household department.

## IMPORTANT
because Julia 1.0 has been published, I am refactoring & transferring this program to Julia 1.0 version. The new version has a better readability and is easy to generalize. New features, such as NamedTuple type and LaTex documentation, will be fully supported in the forthcoming version. And I will udpate the mark system (Greece characters now are widely used as variable names). Readers may find it easier to read, run and generalize the model. However, because the reform of data structure system & that Julia 1.0 is no more based on MKL, the performance of the new version may be affected a little.

by Tianhao Zhao
2018-9-10

(wrote the README in so hurry, pls do not mind my grammar)
