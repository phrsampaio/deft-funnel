# Derivative-free Trust FUNNEL (DEFT-FUNNEL)

This solver searches for the global optimum of the following problem:

min f(x)  
subject to:   
ls <= c(x) <= us,  
lx <=   x  <= ux,  
where f(x) and c(x)=(c_1(x), ..., c_m(x)) are considered as black boxes.

It builds local (at most fully quadratic) interpolation models
from known function values for the objective and constraint functions.
It solves the nonlinear problem using a SQP trust-region-based algorithm
where an active-set method is applied for handling the bound constraints 
and where the convergence is driven by a funnel bound on the constraint 
violation.

In order to find the global optimum, it makes use of a clustering-based 
multi-start technique named Multi-Level Single Linkage (MLSL) to select the 
starting points of the local searches done by the SQP algorithm.

## Author and maintainer: 

Phillipe R. Sampaio  
Veolia Research and Innovation (VERI)  
sampaio.phillipe at gmail.com

## Main references:

* [Ph. L. Sampaio and Ph. L. Toint, 
"Numerical experience with a derivative-free trust-funnel method for nonlinear optimization problems with general nonlinear constraints", Optimization Methods and Software, 31(3), pages 511-534, 2016.](https://www.tandfonline.com/doi/abs/10.1080/10556788.2015.1135919)

* [Ph. L. Sampaio and Ph. L. Toint, 
"A derivative-free trust-funnel method for equality-constrained nonlinear optimization", Computational Optimization and Applications, 61(1), pages 25-49, 2015.](https://doi.org/10.1007/s10589-014-9715-3)

## Contributors: 

Philippe L. Toint (UNamur), Serge Gratton (CERFACS) and 
Anke Troeltzsch (German Aerospace Center, DLR).

## License

This software is released under the MIT license. 
See `License.txt` for more info.

# Using DEFT-FUNNEL without multi-start

Run the 'startup' function to add the necessary folders to the path:
```
>> startup
```

Then call DEFT-FUNNEL at the Matlab command window by typing:
```
>> [ x, fx, norm_c_s, mu, nfeval ] = deft_funnel( @obj_function, @constraints, x0, nbcons )
```

**Mandatory input:**

* *obj_function* : a function handle pointing to the function to be minimized

* *constraints*  : a function handle pointing to contraint functions

* *x0*           : initial guess (no need to be feasible)

* *nbcons*       : number of constraints (bound constraints not included)

**Optional input:**

* *lsbounds*     : vector of lower bounds for the constraints

* *usbounds*     : vector of upper bounds for the constraints

* *lxbounds*     : vector of lower bounds for the x variables

* *uxbounds*     : vector of upper bounds for the x variables

* *maxeval*      : maximum number of evaluations (default: 500*n)

**More parameters:** see deft_funnel_set_parameters.m

**Output:**

* *x*            : the best approximation found to a local minimizer,

* *fx*           : the value of the objective function at x,

* *mu*           : local estimates for the Lagrange multipliers

* *indicators*   : feasibility and optimiality indicators

* *evaluations*  : number of calls to the objective function and constraints

* *iterate*      : info related to the best point found as well as
               the coordinates of all past iterates
               
* *exit_algo*    : output signal (0: terminated with success; -1: terminated with errors)

## Examples of usage:

Some of the test problems included in this package are:

* Problem HS21 (see files problem_hs21_obj.m and problem_hs21_cons.m):
```
>> [ x, fx, mu, indicators, evaluations, iterate, exit_algo ] =      ...
deft_funnel( @problem_hs21_obj, @problem_hs21_cons, [-1 -1], 1,      ...
'lsbounds', 0, 'usbounds', Inf, ...
'lxbounds', [2 -50], 'uxbounds', [50 50] )
```

* Problem HS23 (see files problem_hs23_obj.m and problem_hs23_cons.m):
```
>> [ x, fx, mu, indicators, evaluations, iterate, exit_algo ] =      ...
deft_funnel( @problem_hs23_obj, @problem_hs23_cons, [3 1], 5,        ...
'lsbounds', [0 0 0 0 0], 'usbounds', [Inf Inf Inf Inf Inf],          ...
'lxbounds', [-50 -50], 'uxbounds', [50 50] )
```

A collection of tests may be run through run_deft_funnel.m by typing
```
>> run_deft_funnel
```

at the Matlab command window. run_deft_funnel.m makes use of 3 files that
describe each of the test problems:
1. deft_funnel_problem_init.m (entry parameters), 
2. deft_funnel_problem_obj.m (objective function) 
3. deft_funnel_problem_cons.m (constraint function(s)).

The file deft_funnel_problem_init.m also defines if a constraint in 
deft_funnel_problem_cons.m is an equality or an inequality through 
the lower bounds 'ls' and the upper bounds 'us'.

# Using DEFT-FUNNEL with multi-start

Run the 'startup' function to add the necessary folders to the path:
```
>> startup
```

Then call DEFT-FUNNEL at the Matlab command window by typing:
```
>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @obj_function, @constraints, n, nbcons )
```

**Mandatory input:**

* *obj_function* : a function handle pointing to the function to be minimized

* *constraints*  : a function handle pointing to contraint functions

* *n*            : number of decision variables

* *nbcons*       : number of constraints (bound constraints not included)

No starting point is required from the user in the multi-start case.

**Optional input:**

* *lsbounds*          : vector of lower bounds for the constraints

* *usbounds*          : vector of upper bounds for the constraints

* *lxbounds*          : vector of lower bounds for the x variables

* *uxbounds*          : vector of upper bounds for the x variables

* *maxeval*           : maximum number of evaluations (default: 5000*n)

* *maxeval_ls*        : maximum number of evaluations per local search (default: maxeval*0.7)

* *f_global_optimum*  : known objective function value of the global optimum

**Output:**

* *best_sol*          : best feasible solution found

* *best_fval*         : objective function value of ´best_sol´

* *best_indicators*   : indicators of ´best_sol´

* *total_eval*        : number of evaluations used

* *nb_local_searches* : number of local searches done

* *fX*                : objctive function values of all local minima found

## Examples of usage:

The following problems are included in this package:
```
>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_handbook_quadcons_pb3_obj,                      ...
@problem_handbook_quadcons_pb3_cons, 6, 5,                                       ...
'lsbounds', [ 4 4 -Inf -Inf 2], 'usbounds', [ Inf Inf 2 2 6 ],                   ...
'lxbounds', [ 0 0 1 0 1 0 ], 'uxbounds', [ Inf Inf 5 6 5 10 ])

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_gomez_pb3_obj, @problem_gomez_pb3_cons, 2, 1,   ...
'lsbounds', -Inf, 'usbounds', 0,                                                 ...
'lxbounds', [ -1 -1 ], 'uxbounds', [ 1 1 ])

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_G3_obj, @problem_G3_cons, 2, 1,                 ...
'lsbounds', 0, 'usbounds', 0,                                                    ...
'lxbounds', [ 0 0 ], 'uxbounds', [ 1 1 ])

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_G6_obj, @problem_G6_cons, 2, 2,                 ...
'lsbounds', [ 0 0 ], 'usbounds', [ Inf Inf ],                                    ...
'lxbounds', [ 13 0 ], 'uxbounds', [ 100 100 ])

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_G8_obj, @problem_G8_cons, 2, 2,                 ...
'lsbounds', [ -Inf -Inf ], 'usbounds', [ 0 0 ],                                  ...
'lxbounds', [ 0 0 ], 'uxbounds', [ 10 10 ])

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_G9_obj, @problem_G9_cons, 7, 4,                 ...
'lsbounds', [ -Inf -Inf -Inf -Inf ], 'usbounds', [ 0 0 0 0 ],                    ...
'lxbounds', -10*ones(1,7), 'uxbounds', 10*ones(1,7))

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_G11_obj, @problem_G11_cons, 2, 1,               ...
'lsbounds', 0, 'usbounds', 0,                                                    ...
'lxbounds', -1*ones(1,2), 'uxbounds', ones(1,2))

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_PrW_obj, @problem_PrW_cons, 4, 6,               ...
'lsbounds', [ -Inf -Inf -Inf -Inf -Inf -Inf ], 'usbounds', [ 0 0 0 0 0 0 ],      ...
'lxbounds', [ 0.125 0.1 0.1 0.1 ], 'uxbounds', 10*ones(1,4))

>> [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches, fX ] = ...
deft_funnel_multistart( @problem_PrP_obj, @problem_PrP_cons, 4, 3,               ...
'lsbounds', [ -Inf -Inf -Inf ], 'usbounds', [ 0 0 0 ],                           ...
'lxbounds', [ 0 0 0 0 ], 'uxbounds', [ 1 1 50 240 ] )
```

**CONDITIONS OF USE:** Use at your own risk! No guarantee of any kind given.
