# General Splines Library
Library to represent and formulate motion and trajectory planning problmes with generalized splines.

- Piecewise polynomial curves representation
- Automatic **exact** (and fast) differentiation
- Algebraic operations: inner product, norms, addition, multiplication, composition and concatenation of curves (allows only when it has mathematical sense).
- Optimization with waypoint (via-point) constraints: minimum jerk, snap, crank, etc.

# Motivation
- **Definition** *generalized spline*:
    1. Let <img src="https://render.githubusercontent.com/render/math?math=J=[0, T]"> and consider the partition of  <img src="https://render.githubusercontent.com/render/math?math=J"> given by  <img src="https://render.githubusercontent.com/render/math?math=N+1"> points <img src="https://render.githubusercontent.com/render/math?math=t_i\in J">, i.e. <img src="https://render.githubusercontent.com/render/math?math=I_1, I_2, ... ,I_N"> with <img src="https://render.githubusercontent.com/render/math?math=I_i=[t_i, t_{i+1})">.
    2. Let <img src="https://render.githubusercontent.com/render/math?math=I_0=[-1,1]"> and <img src="https://render.githubusercontent.com/render/math?math=f_1, f_2, ... ,f_k"> be <img src="https://render.githubusercontent.com/render/math?math=k"> linearly independent functions <img src="https://render.githubusercontent.com/render/math?math=f_i:I_0\longrightarrow \mathbb{R}">. 
    3. a

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=s_i(t)= 2\frac{t-t_i}{t_{i+1}-t_i} - 1">
</p>

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

## Implemented basis
Up to now we have implemented tree types of basis
- Basis for the minimum jerk problems, called `cBasis0010`, because their optimize the L2-norm of the jerk
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\Large I=\int_0^T  \left\|\frac{\mathsf{d}^3\mathbf{q}}{\mathsf{d} t^3 }\right\|^2 d t">
</p>

- Basis for the weighed speed-jerk problems, called `cBasis1010`, because their optimize convex combination of the L2-norm of the speed and the L2 norm of the jerk
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\Large I=\int_0^T \alpha  \left\|\frac{\mathsf{d}\mathbf{q}}{\mathsf{d} t }\right\|^2 %2B (\alpha-1)\left\|\frac{\mathsf{d}^3\mathbf{q}}{\mathsf{d} t^3 }\right\|^2 d t">
</p>

# Requirements

- numpy
- scipy
- matplotlib
