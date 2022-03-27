# General Splines Library
Library to represent and formulate motion and trajectory planning problems with generalized splines and piece-wise polynomials.

- Piecewise polynomial curves representation
- Automatic **exact** (and fast) differentiation of generalized splines
- Algebraic operations: inner product, norms, addition, multiplication, composition and concatenation of curves (allows only when it has mathematical sense).
- Optimization with waypoint (via-point) constraints: minimum jerk, snap, crank, etc.
- Implements piece-wise Lagrange polynomials at Gauss-Lobatto points.

# Definition
- **Definition** A **generalized spline** is a piece-wise defined curve such that in each interval it is the linear combination of certain linearly independent functions <img src="https://render.githubusercontent.com/render/math?math=B_1, B_2, ... ,B_k">
- **Formal Definition**
    1. Let <img src="https://render.githubusercontent.com/render/math?math=J=[0, T]"> and consider the partition of  <img src="https://render.githubusercontent.com/render/math?math=J"> given by  <img src="https://render.githubusercontent.com/render/math?math=N %2B 1"> points <img src="https://render.githubusercontent.com/render/math?math=t_i\in J">, i.e. <img src="https://render.githubusercontent.com/render/math?math=I_1, I_2, ... ,I_N"> with <img src="https://render.githubusercontent.com/render/math?math=I_i=[t_i, t_{i %2B 1})">.
    2. Let <img src="https://render.githubusercontent.com/render/math?math=I_0=[-1,1]"> and <img src="https://render.githubusercontent.com/render/math?math=B_1, B_2, ... ,B_k"> be <img src="https://render.githubusercontent.com/render/math?math=k"> linearly independent functions <img src="https://render.githubusercontent.com/render/math?math=B_i:I_0\longrightarrow \mathbb{R}">. 
    3. Let <img src="https://render.githubusercontent.com/render/math?math=s_i:I_i\longrightarrow I_0"> given by
    <p align="center">
    <img src="https://render.githubusercontent.com/render/math?math=s_i(t)= 2\frac{t-t_i}{t_{i %2B 1}-t_i} - 1">
    </p>

    4. Let <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}_i^j \in\mathbb{R}^n">.

    5. A **generalized spline** from <img src="https://render.githubusercontent.com/render/math?math=J"> into <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^n"> is a curve given by
    <p align="center">
    <img src="https://render.githubusercontent.com/render/math?math=f=\sum_{i=1}^k \mathbf{y}_k^j B_k \circ s_j(t) \text{ if } t\in I_j\ \ \ \ \ \ \ \ \ \ \ \ \ (\star)">
    </p>

# Motivation
Generalized splines appear naturally in problems of trajectory optimization when waypoint constraints are added.
In other words, if we desire to optimize a motion which pass trough a sequence of positions we will meet with generalized splines.

Generalized splines trajectories arising in such kind of optimization problems are uniquely characterized by the waypoints that they attain, the time intervals between waypoints and the speed, and possible higher order derivatives at the boundaries.
Moreover, such a relation is synthesised in the expression of the type

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\mathbf{A}(\boldsymbol\tau)\mathbf{y} = \mathbf{b}(\mathbf{w})\ \ \ \ \ \ \ \ \ \ \ \ \ \ (0)">
</p>

where <img src="https://render.githubusercontent.com/render/math?math=\mathbf{A}(\boldsymbol\tau)"> is a matrix which depends on the time intervals <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol\tau">, <img src="https://render.githubusercontent.com/render/math?math=\mathbf{b}(\mathbf{w})"> is a column vector which depends on the waypoints <img src="https://render.githubusercontent.com/render/math?math=\mathbf{w}">, the speed, and possible higher order derivatives at the boundaries, and <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}"> is a column vector which represents uniquely the curve at each interval.

The main challenge to build into a computer a trajectory optimization problems with waypoint constraints is to compute the derivatives of <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}"> with respect to <img src="https://render.githubusercontent.com/render/math?math=\mathbf{\tau}">. 

This library provides a uniform and simple interface to formulate gradient-based optimization problems for waypoint-constrained trajectory planning. The library leverage on the representation (0) to compute the "derivatives of the splines" with respect to the time intervals (and possibly the waypoints) as the corresponding derivatives of <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}">


# Background

This library is aimed to find a trajectory passing trough a sequence of waypoints <img src="https://render.githubusercontent.com/render/math?math=\{\mathbf{w}_0, ...,\mathbf{w}_{N %2B 1}\}"> such that the following integral is minized
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\Large I=\int_0^T \alpha_1\left\|\frac{\mathsf{d}\mathbf{q}}{\mathsf{d} t }\right\|^2 %2B \alpha_2 \left\|\frac{\mathsf{d}^2\mathbf{q}}{\mathsf{d} t^2 }\right\|^2 %2B \alpha_3\left\|\frac{\mathsf{d}^3\mathbf{q}}{\mathsf{d} t^3 }\right\|^2 %2B  \alpha_4\left\|\frac{\mathsf{d}^4\mathbf{q}}{\mathsf{d} t^4 }\right\|^2 \mathsf{d} t \ \ \ \ \ (1)">
</p>
It may be proven that such a problem can be subdivided in two steps

 1. Find the family of optimal curves that joint waypoints
 2. Compute time instants <img src="https://render.githubusercontent.com/render/math?math=\{t_0, t_1,  ...,t_N, t_{N %2B 1}\}"> where the optimal curves must be attached

The step 1. is done by solving a linear ordinary differential equation. One method to achieve 2. is to formulate an optimization problem (e.g. a gradient based one).

## Optimal curves
We underline that this library leverages on the [general theory of linear ODEs](https://en.wikipedia.org/wiki/Linear_differential_equation).
It may be proven that any optimal of (1) solves the following linear ODE, which turn out to be the Euler-Lagrange equations at each interval <img src="https://render.githubusercontent.com/render/math?math=[t_i, t_{i %2B 1}]">
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=-\alpha_1\frac{\mathsf{d}^2\mathbf{q}}{\mathsf{d} t^2 } %2B \alpha_2 \frac{\mathsf{d}^4\mathbf{q}}{\mathsf{d} t^4 } - \alpha_3\frac{\mathsf{d}^6\mathbf{q}}{\mathsf{d} t^6 } %2B  \alpha_4 \frac{\mathsf{d}^8\mathbf{q}}{\mathsf{d} t^8 } = 0\ \ \ \ \ (2)">
</p>
with the following boundary conditions
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\mathbf{q}(t_i) = \mathbf{w}_i\ \ \ \ \ \ \ \mathbf{q}(t_{i%2B1}) = \mathbf{w}_{i%2B1}\ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)">
</p>
Because the ODE (2) is linear, we can compute its general suction depending on the value of the coefficients <img src="https://render.githubusercontent.com/render/math?math=\alpha_i">.

In fact, the general solution of (2) may be written as a piecewise function defined at each interval as
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\mathbf{q} = \sum_{i=1}^{n_b} \mathbf{y}_i^j B_i(t) \ \ \ \ \text{if}\ \ \ \ t \in [t_{j}, t_{j %2B 1}]\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (4)">
</p>

where <img src="https://render.githubusercontent.com/render/math?math=n_b">, <img src="https://render.githubusercontent.com/render/math?math=B_i(t)"> depend on the coefficients <img src="https://render.githubusercontent.com/render/math?math=\alpha_i"> and <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}_i^j"> are column vectors in which represents the curve uniquely at the interval <img src="https://render.githubusercontent.com/render/math?math=[t_j, t_{j %2B 1}]">.

If we stack the column vectors <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}_i^j"> in a suitable way we obtain the column vector <img src="https://render.githubusercontent.com/render/math?math=\mathbf{y}"> used in (0). In fact, the equation (0) is obtained after applying to (4) the waypoint constrains and the boundary conditions.

After substituting (4) in (1) we obtain the following expression
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=I=\mathbf{y}^\top \mathbf{Q}(\boldsymbol{\tau}) \mathbf{y}\ \ \ \ \ \ \ \ \ \ \ \ \ (5)">
</p>

Finally we can substitute (0) in (5) to obtain the representation of (1) subject to the waypoint constrains as a function of <img src="https://render.githubusercontent.com/render/math?math=N"> real variables:

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=I=I(\boldsymbol{\tau})=\mathbf{b}^{\top}\mathbf{A}^{-\top}(\boldsymbol{\tau})\mathbf{Q}(\boldsymbol{\tau}) \mathbf{A}^{-1}(\boldsymbol{\tau})\mathbf{b}\ \ \ \ \ \ \ \ \ \ \ \ \ (6)">
</p>

# Software architecture

From the formalization of the optimization problem, we derive that a flexible and uniform methodology for the construction of the problem of optimizing (1) consists in designing an abstract representation of the basis <img src="https://render.githubusercontent.com/render/math?math=B_i(t)"> in (4) capable of building in an automatic fashion the constraint (0), the cost function (5) and their derivatives.

In fact, note that the input of any gradient-based optimizer is the expression (6)  and its derivatives. 
This library provides a template class to represent the basis <img src="https://render.githubusercontent.com/render/math?math=B_i(t)"> and a series of procedures which utilizes these basis as an input and then generate (0), (6) and their derivatives.

This library provides the class `gsplines::basis::Basis` which represent an arbitrary set of linearly independent functions and the class `gsplines::GSpline` which implements the definition <img src="https://render.githubusercontent.com/render/math?math=\star">.
In addition this library provides the class `gsplines::functions::FunctionBase` and `gsplines::functions::FunctionExpression` that allows to define arbitrary curves in <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^n"> with the possibility of performing algebraic operations and derivation.

## Arbitrary function definition and algebraic operations

The abstract class [`gsplines::functions::FunctionBase`](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/master/include/gsplines/Functions/FunctionBase.hpp) represents objects that can be evaluated at points of some interval of <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}">. To inherit from this class the auxiliary template class [`gsplines::functions::FunctionInheritanceHelper`](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/master/include/gsplines/Functions/FunctionInheritanceHelper.hpp) allows to define a custom function which <img src="https://render.githubusercontent.com/render/math?math=k">-derivative is implemented as another class using the [Curiously recurring template pattern](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern). For example, the exponential and other elementary functions are [declared here](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/fa7395af1719a9c560eba15ff4bc68584da2c7e7/include/gsplines/Functions/ElementalFunctions.hpp#L82).


The class `gsplines::functions::FunctionExpression` inherits from `gsplines::functions::FunctionBase` and contain an array of pointers to `gsplines::functions::FunctionBase` called `function_array_`. The evaluation operator of `gsplines::functions::FunctionExpression` evaluate the desired algebraic operation between the functions in `function_array_`.
This architecture allow to implement complex algebraic expressions by recursively calling the evaluation of each function in `function_array_`.


# Requirements

- `numpy`
- `scipy`
- `matplotlib`
- `coinor-libipopt-dev`
- `if-opt`

# Publications

- Rojas, Rafael A., et al. "A variational approach to minimum-jerk trajectories for psychological safety in collaborative assembly stations." IEEE Robotics and Automation Letters 4.2 (2019): 823-829.
```
@article{rojas2019variational,
  title={A variational approach to minimum-jerk trajectories for psychological safety in collaborative assembly stations},
  author={Rojas, Rafael A and Garcia, Manuel A Ruiz and Wehrle, Erich and Vidoni, Renato},
  journal={IEEE Robotics and Automation Letters},
  volume={4},
  number={2},
  pages={823--829},
  year={2019},
  publisher={IEEE}
}
```

- Rojas, Rafael A., et al. "Combining safety and speed in collaborative assembly systems–An approach to time optimal trajectories for collaborative robots." Procedia CIRP 97 (2021): 308-312.
```
@article{rojas2021combining,
  title={Combining safety and speed in collaborative assembly systems--An approach to time optimal trajectories for collaborative robots},
  author={Rojas, Rafael A and Garcia, Manuel A Ruiz and Gualtieri, Luca and Rauch, Erwin},
  journal={Procedia CIRP},
  volume={97},
  pages={308--312},
  year={2021},
  publisher={Elsevier}
}
```

- Rojas, Rafael A., and Renato Vidoni. "Designing fast and smooth trajectories in collaborative workstations." IEEE Robotics and Automation Letters 6.2 (2021): 1700-1706.
```
@article{rojas2021designing,
  title={Designing fast and smooth trajectories in collaborative workstations},
  author={Rojas, Rafael A and Vidoni, Renato},
  journal={IEEE Robotics and Automation Letters},
  volume={6},
  number={2},
  pages={1700--1706},
  year={2021},
  publisher={IEEE}
}
```

- Rojas, Rafael A., Erich Wehrle, and Renato Vidoni. "A multicriteria motion planning approach for combining smoothness and speed in collaborative assembly systems." Applied Sciences 10.15 (2020): 5086.
```
@article{rojas2020multicriteria,
  title={A multicriteria motion planning approach for combining smoothness and speed in collaborative assembly systems},
  author={Rojas, Rafael A and Wehrle, Erich and Vidoni, Renato},
  journal={Applied Sciences},
  volume={10},
  number={15},
  pages={5086},
  year={2020},
  publisher={Multidisciplinary Digital Publishing Institute}
}
```
