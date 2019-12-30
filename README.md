# Derivative-free Trust FUNNEL

This solver searches for the global minima of grey-box and black-box 
optimization problems as defined below:

min f(x)  
subject to:   
lc <= c(x) <= uc,  
lh <= h(x) <= uh,  
lx <=   x  <= ux,  
where f(x) might be or not a black box, c(x)=(c_1(x), ..., c_q(x)) are black-box 
constraint functions, h(x)=(h_1(x), ..., h_l(x)) are white-box constraint functions 
(i.e. their analytical expressions as well as their derivatives are available), 
lc, lh, uc and uh are vectors defining the lower and upper bounds of c(x) and h(x), 
and lx and ux are lower and upper bounds on x.

DEFT-FUNNEL builds local (at most fully quadratic) interpolation models
from known function values for the black-box functions. It solves the 
nonlinear problem using a SQP trust-region-based algorithm where an active-set 
method is applied for handling the bound constraints and where the convergence 
is driven by a funnel bound on the constraint violation.

In order to find the global minimum, it makes use of a clustering-based 
multistart technique called Multi-Level Single Linkage (MLSL) to select the 
starting points of the local searches done by the SQP algorithm.

## Author and maintainer 

Phillipe Rodrigues Sampaio  
Veolia Research and Innovation (VERI)  
sampaio.phillipe at gmail.com

## Main references

* [Ph. R. Sampaio and Ph. L. Toint, 
"Numerical experience with a derivative-free trust-funnel method for nonlinear optimization problems with general nonlinear constraints", Optimization Methods and Software, 31(3), pages 511-534, 2016.](https://www.tandfonline.com/doi/abs/10.1080/10556788.2015.1135919)

* [Ph. R. Sampaio and Ph. L. Toint, 
"A derivative-free trust-funnel method for equality-constrained nonlinear optimization", Computational Optimization and Applications, 61(1), pages 25-49, 2015.](https://doi.org/10.1007/s10589-014-9715-3)

## Contributors 

Philippe L. Toint (UNamur), Serge Gratton (CERFACS) and 
Anke Troeltzsch (German Aerospace Center, DLR).

## License

This software is released under the MIT license. 
See `LICENSE.md` for more info.

# DEFT-FUNNEL without multistart

If global minima are not required and local minima are enough, call DEFT-FUNNEL 
at the Matlab command window by typing:
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] = deft_funnel(@f, @c, @h, ...
@dev_f, @dev_h, x0, nb_cons_c, nb_cons_h)
```

**Mandatory input:**

* *f*            : function handle of the objective function

* *c*            : function handle of the black-box constraints if any or an empy array [] otherwise

* *h*            : function handle of the white-box constraints if any or an empy array [] otherwise

* *dev_f*        : function handle of the derivatives of f if it is a white box or an empy array [] otherwise

* *dev_h*        : function handle of the derivatives of h if any or an empy array [] otherwise

* *x0*           : starting point (no need to be feasible)

* *nb_cons_c*    : number of black-box constraints (bound constraints not included)

* *nb_cons_h*    : number of white-box constraints (bound constraints not included)

**IMPORTANT**:

The output of 'dev_f' must be a cell array containing two cells, one for each component below: 

* gf             : the gradient of 'f' with dimensions 'n' x 1.

* Hf             : the hessian matrix with dimensions 'n' x 'n'.

The output of 'dev_h' must be a cell array containing two cells, one for each component below: 

* Jh             : the Jacobian matrix of 'h' with dimensions 'nb_cons_h' x 'n'.

* Hh             : a matrix containing the hessians of each constraint put 
side by side, i.e., a matrix with dimensions 'n' x ('nb_cons_h' * 'n').

**Optional input:**

* *lsbounds*     : vector of lower bounds for all the constraints (white box and black box)

* *usbounds*     : vector of upper bounds for all the constraints (white box and black box)

* *lxbounds*     : vector of lower bounds for the x variables

* *uxbounds*     : vector of upper bounds for the x variables

* *maxeval*      : maximum number of evaluations (default: 500*n)

* *type_f*       : string 'BB' if f is a black box (default) or 'WB' otherwise

* *whichmodel*   : approach to build the surrogate models: 

    - 0: Subbasis model
    
    - 1: Frobenius-norm model
    
    - 2: minimum l2-norm (Default)
    
    - 3: regression (recommended for noisy functions)

**More parameters:** see 'deft_funnel_set_parameters.m'.

**Output:**

* *x*            : the best approximation found to a local minimum,

* *fx*           : the value of the objective function at x,

* *mu*           : local estimates for the Lagrange multipliers

* *indicators*   : feasibility and optimality indicators

* *evaluations*  : number of calls to the objective function and constraints

* *iterate*      : more info related to the best point found as well as
                   the coordinates of all past iterates
               
* *exit_algo*    : output signal (0: terminated with success; -1: terminated with errors)

## Examples of usage

If 'f' is a black box, 'dev_f' is expected to be an empty array:

```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] = deft_funnel(@f, @c, @h, ...
[], @dev_h, x0, nb_cons_c, nb_cons_h)
```

If 'f' is a white box, the user must indicate it (type_f='WB') through the input 
argument 'type_f' as in the example below (by default, type_f='BB'). 
In this case, 'dev_f' is expected to be a valid function that computes 
the derivatives of 'f'.

```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] = deft_funnel(@f, @c, @h, ...
@dev_f, @dev_h, x0, nb_cons_c, nb_cons_h, 'type_f', 'WB')
```

If there are no black-box constraints but only white-box ones, an empty array 
must be used in the place of @c and 'nb_cons_c' must equal 0:

```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] = deft_funnel(@f, [], @h, ...
@dev_f, @dev_h, x0, 0, nb_cons_h, 'type_f', 'WB')
```

If there are no white-box constraints but only black-box ones, an empty array 
must be used in the place of @h and @dev_h; moreover, 'nb_cons_h' must equal 0:

```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] = deft_funnel(@f, @c, [], ...
@dev_f, [], x0, nb_cons_c, 0)
```

## Black-box test problems

Some of the black-box (BB) test problems included in this package are:

* Problem HS21 (see files problem_hs21_obj.m and problem_hs21_cons.m):
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
deft_funnel(@problem_hs21_obj, @problem_hs21_cons,             ...
[], [], [], [-1 -1], 1, 0, 'lsbounds', 0, 'usbounds', Inf,     ...
'lxbounds', [2 -50], 'uxbounds', [50 50])
```

* Problem HS23 (see files problem_hs23_obj.m and problem_hs23_cons.m):
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
deft_funnel(@problem_hs23_obj, @problem_hs23_cons,             ...
[], [], [], [3 1], 5, 0, 'lsbounds', [0 0 0 0 0],              ...
'usbounds', [Inf Inf Inf Inf Inf],                             ...
'lxbounds', [-50 -50], 'uxbounds', [50 50])
```

A collection of 9 BB test problems may be solved by typing
```
>> startup
```
within the 'testset' directory and then
```
>> run_deft_funnel_all_bb_test_probs
```
All the BB test problems are found in the directory 'testset/blackbox'. 
The 'startup' function adds the folder of test problems to the path and 
'run_deft_funnel_all_bb_test_probs' calls DEFT-FUNNEL within a loop to solve those 
problems.

In order to solve a specific test problem, type
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
run_deft_funnel_single_bb_test_prob(nprob)
```
where 'nprob' is a number from 1 to 9.

The files 'run_deft_funnel_all_bb_test_probs.m' and 
'run_deft_funnel_single_bb_test_prob.m' make use of 6 other files that
describe the test problems:
1. deft_funnel_problem_init.m - entry parameters, 
2. deft_funnel_problem_obj.m - objective function,
3. deft_funnel_problem_cons_c.m - black-box constraint function(s).
4. deft_funnel_problem_cons_h.m - white-box constraint function(s).
5. deft_funnel_problem_dev_f.m - derivatives of the objective function.
6. deft_funnel_problem_dev_h.m - derivatives of the white-box constraints.

The file 'deft_funnel_problem_init.m' also defines if a constraint in 
'deft_funnel_problem_cons_c.m' and 'deft_funnel_problem_cons_h.m' is an equality or 
an inequality through the lower bounds 'ls' and the upper bounds 'us'. 

Finally, the test problems defined in 'deft_funnel_problem_init.m' ranging 
from 10 to 28 are designed for global optimization, so they should be solved 
with multistart. In particular, the problems 24-28 are of grey-box type.

# DEFT-FUNNEL with multistart

Running DEFT-FUNNEL with multistart is recommended in the following scenarii:
* global minima are required;
* the objective function is known to be multimodal;
* previous trials with a specific starting point were not successful;
* the shape of the objective function is unknown but something in your heart says that it is nonlinear.

DEFT-FUNNEL with multistart can be called by typing:
```
>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@f, @c, @h, @dev_f, @dev_h, n, nb_cons_c, nb_cons_h)
```

**Mandatory input:**

* *f*            : function handle of the objective function

* *c*            : function handle of the black-box constraints if any or an empy array [] otherwise

* *h*            : function handle of the white-box constraints if any or an empy array [] otherwise

* *dev_f*        : function handle of the derivatives of f if it is a white box or an empy array [] otherwise

* *dev_h*        : function handle of the derivatives of h if any or an empy array [] otherwise

* *n*            : number of decision variables

* *nb_cons_c*    : number of black-box constraints (bound constraints not included)

* *nb_cons_h*    : number of white-box constraints (bound constraints not included)

No starting point is required from the user in the multistart case.

**IMPORTANT**:

The output of 'dev_f' must be a cell array containing two cells, one for each component below: 

* gf             : the gradient of 'f' with dimensions 'n' x 1.

* Hf             : the hessian matrix with dimensions 'n' x 'n'.

The output of 'dev_h' must be a cell array containing two cells, one for each component below: 

* Jh             : the Jacobian matrix of 'h' with dimensions 'nb_cons_h' x 'n'.

* Hh             : a matrix containing the hessians of each constraint put 
side by side, i.e., a matrix with dimensions 'n' x ('nb_cons_h' * 'n').

**Optional input:**

* *lsbounds*          : vector of lower bounds for all the constraints (white box and black box)

* *usbounds*          : vector of upper bounds for all the constraints (white box and black box)

* *lxbounds*          : vector of lower bounds for the x variables

* *uxbounds*          : vector of upper bounds for the x variables

* *maxeval*           : maximum number of evaluations (default: 5000*n)

* *maxeval_ls*        : maximum number of evaluations per local search (default: maxeval*0.7)

* *type_f*            : string 'BB' if f is a black box (default) or 'WB' otherwise

* *whichmodel*        : approach to build the surrogate models: 

    - 0: Subbasis model
    
    - 1: Frobenius-norm model
    
    - 2: minimum l2-norm (Default)
    
    - 3: regression (recommended for noisy functions)

* *f_global_optimum*  : known objective function value of the global optimum

**Output:**

* *best_sol*          : best feasible solution found

* *best_feval*        : objective function value of ´best_sol´

* *best_indicators*   : indicators of ´best_sol´

* *total_eval*        : number of evaluations used

* *nb_local_searches* : number of local searches done

* *fL*                : objective function values of all local minima found

## Examples of usage

See the examples for the case without multistart as the only difference 
is that the user must pass the number of decision variables rather
than the starting point as input.

## Black-box test problems

A collection of well-known test problems for constrained global optimization 
are included in this package. In the examples below, all the functions are 
black boxes. In order to run DEFT-FUNNEL on any of those problems, 
first type:
```
>> startup
```
within the 'testset' directory. 

Some of the test problems collected from the literature are the following:
```
>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_hesse_obj,                                      ...
@problem_hesse_cons, [], [], [], 6, 5, 0,                                       ...
'lsbounds', [4 4 -Inf -Inf 2], 'usbounds', [Inf Inf 2 2 6],                     ...
'lxbounds', [0 0 1 0 1 0], 'uxbounds', [Inf Inf 5 6 5 10])

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_gomez_pb3_obj, @problem_gomez_pb3_cons,         ...
[], [], [], 2, 1, 0, 'lsbounds', -Inf, 'usbounds', 0,                           ...
'lxbounds', [-1 -1], 'uxbounds', [1 1])

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_G3_obj, @problem_G3_cons,                       ...
[], [], [], 2, 1, 0, 'lsbounds', 0, 'usbounds', 0,                              ...
'lxbounds', [0 0], 'uxbounds', [1 1])

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_G6_obj, @problem_G6_cons,                       ...
[], [], [], 2, 2, 0, 'lsbounds', [0 0], 'usbounds', [Inf Inf],                  ...
'lxbounds', [13 0], 'uxbounds', [100 100])

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_G8_obj, @problem_G8_cons,                       ...
[], [], [], 2, 2, 0, 'lsbounds', [-Inf -Inf], 'usbounds', [0 0],                ...
'lxbounds', [0 0], 'uxbounds', [10 10])

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_G9_obj, @problem_G9_cons,                       ...
[], [], [], 7, 4, 0, 'lsbounds', [-Inf -Inf -Inf -Inf],                         ...
'usbounds', [0 0 0 0], 'lxbounds', -10*ones(1,7),                               ...
'uxbounds', 10*ones(1,7))

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_G11_obj, @problem_G11_cons,                     ...
[], [], [], 2, 1, 0, 'lsbounds', 0, 'usbounds', 0,                              ...
'lxbounds', -1*ones(1,2), 'uxbounds', ones(1,2))

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_WB4_obj, @problem_WB4_cons,                     ...
[], [], [], 4, 6, 0, 'lsbounds', [-Inf -Inf -Inf -Inf -Inf -Inf],               ...
'usbounds', [0 0 0 0 0 0], 'lxbounds', [0.125 0.1 0.1 0.1],                     ...
'uxbounds', 10*ones(1,4))

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@problem_PVD4_obj, @problem_PVD4_cons,                   ...
[], [], [], 4, 3, 0, 'lsbounds', [-Inf -Inf -Inf], 'usbounds', [0 0 0],         ...
'lxbounds', [0 0 0 0], 'uxbounds', [1 1 50 240])
```

The collection of 23 BB test problems may be solved with multistart by typing
```
>> startup
```
within the 'testset' directory and then
```
>> run_deft_funnel_multistart_all_bb_test_probs
```
In order to solve a specific test problem with multistart, type
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
run_deft_funnel_multistart_single_bb_test_prob(nprob)
```
where 'nprob' is a number from 1 to 23.

## Grey-box test problems

The 5 grey-box test problems included here were constructed from 
the BB test problems either by making the objective function white box or 
by transforming some of the BB constraints into white boxes. They are found 
at the directory 'testset/greybox'. 

In order to run DEFT-FUNNEL on them, you must type 
```
>> startup
```
within the 'testset' directory. 

Two of the 5 test problems available are:

* Problem HS21 with the objective function as white box (see files 
'problem_greybox_hs21_obj.m', 'problem_greybox_hs21_cons_c.m' and 
'problem_greybox_hs21_dev_f.m'):
```
[x, fx, mu, indicators, evaluations, iterate, exit_algo] =                  ...
deft_funnel_multistart(@problem_greybox_hs21_obj,                           ...
@problem_greybox_hs21_cons_c, [], @problem_greybox_hs21_dev_f, [],          ...
2, 1, 0, 'lsbounds', 0, 'usbounds', Inf, 'lxbounds', [2 -50],               ...
'uxbounds', [50 50], 'type_f', 'WB')
```

* Problem HS23 with some of the constraints as white boxes (see files 
'problem_greybox_hs23_obj.m', 'problem_greybox_hs23_cons_c.m', 
'problem_greybox_hs23_cons_h.m' and 'problem_greybox_hs23_dev_h.m'):
```
[x, fx, mu, indicators, evaluations, iterate, exit_algo] =                  ...
deft_funnel_multistart(@problem_greybox_hs23_obj,                           ...
@problem_greybox_hs23_cons_c, @problem_greybox_hs23_cons_h,                 ...
@problem_greybox_hs23_dev_f, @problem_greybox_hs23_dev_h, 2, 2, 3,          ...
'lsbounds', [0 0 0 0 0], 'usbounds', [Inf Inf Inf Inf Inf],                 ...
'lxbounds', [-50 -50], 'uxbounds', [50 50], 'type_f', 'WB')
```

In order to solve a specific grey-box test problem with multistart, type
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
run_deft_funnel_multistart_single_gb_test_prob(nprob)
```
where 'nprob' is a number from 24 to 28.

# Evaluation of objective and black-box constraints from a single black-box call

In this case, the first argument must contain the black-box function handle and 
the second argument must be the string 'combined' as in the examples below:
```
>> [x, fx, mu, indicators, evaluations, iterate, exit_algo] =                  ...
deft_funnel(@blackbox, 'combined', @h, [], @dev_h, x0, nb_cons_c, nb_cons_h)

>> [best_sol, best_feval, best_indicators, total_eval, nb_local_searches, fL] = ...
deft_funnel_multistart(@blackbox, 'combined', @h, [], @dev_h, n, nb_cons_c, nb_cons_h)
```

The solver assumes that the first output of @blackbox contains the objective 
function evaluation while the remaining outputs are the BB constraint functions 
evaluations.

**CONDITIONS OF USE:** Use at your own risk! No guarantee of any kind given.
