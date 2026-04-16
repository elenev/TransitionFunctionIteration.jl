# TransitionFunctionIteration.jl

This is a package for globally solving macroeconomic models using Transition Function Iteration of Elenev, Landvoigt, and Van Nieuwerburgh (2021) or, optionally, also using the Simultaneous Transition and Policy Function Iteration (STPFI) algorithm introduced in Cao, Luo, and Nie (2023).

A solved model is represented by three sets of functions of state variables, linearly interpolated on a (potentially irregular) rectangular grid. The first state is always the index of the Markov Chain governing exogenous state variable dynamics.
- Solution functions `Pfct` containing the solutions to a nonlinear system of equations defining the equilibrium at each grid point in the state space.
- Next-period functions `Vfct` containing equilibrium quantities and prices whose next-period values enter in expectation equilibrium conditions today. 
- Transition functions `Tfct` which specify the next-period value of endogenous variables for each realization of exogenous variables next period.

The algorithm initializes these functions based on their deterministic steady-state values and iterates to convergence. The package also provides utilities to simulate, compute impulse responses, plot, and compare results.

## Installation

Install the package using:

```julia-repl
julia> using Pkg
julia> Pkg.add("TransitionFunctionIteration")
```

## Usage 

### Pre-requisites 
To use, you must (1) define a custom struct that inherits from `AbstractModel` and contains your model parameters, and (2) implement a set of functions containing the "economics" of the model. The easiest way to get started is to run:

```julia-repl
julia> generate_template("MyModel","/path/to/model_file.jl")
```

This will create a template file with shells for the functions you need to implement. They are:

- `Exogenv(mdl::MyModel)` to overload the constructor of the exogenous environment struct `Exogenv` using the parameters of the model.
- `define_endog_grids(mdl::MyModel, ::AbstractSteadyState)` to provide grid points for endogenous state variables
- `compute_next_period_values(next, mdl::MyModel, mode::EvaluationMode)::NextPeriod` to provide, for a given exogenous state next period, elements of expectations operators in the equilibrium conditions, values of next-period variables, which enter transition functions, and any other next-period variables, whose values are needed to compute quantities of interest in simulation for a solved model.
- `calc_eqm(state, policy, mdl::MyModel, np::NextPeriod, mode::EvaluationMode)` to evaluate the system of equations, whose solution defines an equilibrium at a given grid point, and to provide updated `Vfct` and `Tfct` values at the solution.

### Usage

Having defined the functions and instantiated your model `mdl = MyModel()`, solve it by:

1. Calculating the steady state and initializing the dynamic model solution at the steady state:

```julia-repl
julia> obj, stv = main_setup()
```

2. Solve the model:

``` julia-repl
julia>main_execute!(obj)
```

Subsequently, you can simulate the model using `simulate`, summarize the simulation using `describe`, and calculate impulse reponses using `simulate_irf`. An array of plotting tools is also available.

### Extending Functionality

Many other functions can be overloaded to extend functionality. Their shells are commented out. See below for details.

- One of, depending on the method of computing the steady state:
    - `steadystate` to define the deterministic steady state analytically (using `AnalyticalSteadyStateMode`)
    - `steadystate_system` to define the deterministic steady state as a solution to a different (usually smaller) system of equations than the dynamic system (using `SystemSteadyStateMode`)
 Alternatively, the steady state can be calculated by solving the dynamic system of equations augmented with constant laws of motion for endogenous state variables. To do this, provide an initial guess when constructing `FixedPointSteadyStateMode`

- `initial_guess_from_stv(stv::SteadyState, mdl::MyModel)` to adjust the steady state values before they're used to initialize the dynamic system.

- `normalize_next_period_values(np::NextPeriod, mdl::MyModel, mode::EvaluationMode)` to transform next-period expectations before they enter the equilibrium conditions in `calc_eqm`. This can be useful if the transformation does not depend on the unknowns being solved for in `calc_eqm`, because it ensures that the (potentially costly) calculations involved in these transformations are done once, and not repeated for every guess as the system is being solved.

- `normalize_fx_errors(::MyModel, fx_nt, simvars)` to normalize errors/residuals of equilibrium conditions, to make them comparable.

- `augment_simulation(sim::StructVector, mdl::MyModel)` to augment a single $T$-period simulation of the model with additional variables not defined in `calc_eqm`. This is especially useful when these variables are _dynamic_, i.e., when they use values from adjacent periods (e.g., returns, growth rates).

## Examples

The examples folder contains two models solved using this package.

1. `ngm.jl`: Neoclassical Growth Model. The textbook model defined by equations below. This is an easy model to solve and, the log-linear approximation is very close to the global solution. In fact, the packages test suite uses the log-linear solution as a "check" on the correct implementation of the method.

$$
\begin{aligned}
 V(A,K) &= \max_{C,K'} \left\{ \frac{C^{1-\gamma}}{1-\gamma} + \beta E[V(A',K') | A] \right\} \\
 &\text{s.t.} \\
 C + I &= A K^\alpha \\
 K' &= (1-\delta) K + I
\end{aligned}
$$

2. `intermedap.jl`: An intermediary-asset pricing model with two agents, occasionally binding constraints, trend growth, and Epstein-Zin preferences. `examples/intermedap.pdf` contains the model details.
