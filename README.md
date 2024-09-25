# General Splines: A Library for Algebraically and Analytically Consistent Representation of Motions
Library to represent and formulate motion and trajectory planning problems with generalized splines and piece-wise polynomials.

- Generalized splines as a `GSpline` class. They can represent
    - [Piecewise](https://en.wikipedia.org/wiki/Piecewise_function) polynomial curves representation
    - [Piecewise](https://en.wikipedia.org/wiki/Piecewise_function) Lagrange polynomials at (interpolation at Gauss-Lobatto points already implemented).
    - **Analytically consistency** `GSpline` provide a `derivate` method which returns its derivative as a new `GSpline` instance. This library provides automatic **exact** ($\matcal{O}(n^2)$) differentiation of the generalized splines implemented.
    - **Algebraic consistency**: This library implements basic operations between `GSplines` inner product, norms, addition, multiplication, composition and concatenation of curves (allows only when it has mathematical sense).
- Optimization with waypoint (via-point) constraints: minimum jerk, snap, crank, etc.
- **ROS 1 implementation** [visit this link](https://github.com/rafaelrojasmiliani/gsplines_cpp_ros)
- **MoveIt implementation** [visit this link](https://github.com/rafaelrojasmiliani/gsplines_moveit)
- **Minimum time path-consistent stop** [visit this link](https://github.com/rafaelrojasmiliani/opstop_cpp)
- Contact: Rafael A. Rojas rafaelrojasmiliani@gmail.com
- **Docker containers with this library installed**
    - *vim awesome plugins for development and moveit*  rafa606/moveit-gsplines-vim-dev:noetic
    - *vim awesome plugins for development and awesome ros packages*  rafa606/ros-gsplines-vim-dev:noetic


# Examples

Get your minimum-X  trajectory passing through random waypoints
```PYTHON
import numpy as np
import gsplines
import gsplines.plot
import gsplines.optimization


numberOfWaypoints = 4
dimensionOfAmbientSpace = 2

waypoints = np.random.rand(numberOfWaypoints, dimensionOfAmbientSpace)

c1 = gsplines.optimization.broken_lines_path(waypoints)
c2 = gsplines.optimization.minimum_acceleration_path(waypoints)
c3 = gsplines.optimization.minimum_jerk_path(waypoints)
c4 = gsplines.optimization.minimum_snap_path(waypoints)
c5 = gsplines.optimization.minimum_crackle_path(waypoints) # Numerically unstable !

gsplines.plot.plot2d_compare([c1, c2, c3, c4, c5], [
                             'green', 'blue', 'magenta', 'red', 'black'],
                             ['min vel', 'min acceleration', 'min jerk',
                              'min snap', 'min crackle'])

```

# Installation

## In Ubuntu using deb packages and ROS
The installation using debian packages requires to have access to the ROS 1 repos ([visit this link](http://wiki.ros.org/it/noetic/Installation/Ubuntu)).
The reason to use ros packages is that this library depends on [`ifopt`](https://github.com/ethz-adrl/ifopt), and its deb package is available with ros.
1. Install the requirements
```bash
sudo apt-get install  python3-matplotlib ros-noetic-ifopt
```
2. Download the package a install
```bash
wget https://github.com/rafaelrojasmiliani/gsplines_cpp/releases/download/master/gsplines-0.0.1-amd64.deb
sudo dpkg -i gsplines-0.0.1-amd64.deb
```

## From source without ros

1. Install the requirements
```bash
sudo apt-get install  python3-matplotlib libgtest-dev cmake libeigen3-dev coinor-libipopt-dev python3-dev
```

2. Install `ifopt`
```bash
   git clone https://github.com/ethz-adrl/ifopt.git
   cd ifopt
   mkdir build
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=/usr
   make -j$(nproc)
   sudo make install
```
3. Download the repo with recursive mode and compile
```bash
git clone --recursive https://github.com/rafaelrojasmiliani/gsplines_cpp.git
cd gsplines_cpp
mkdir build
cd build
cmake .. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/usr
make -j$(nproc)
sudo make install
```

# Definition
- **Definition** A **generalized spline** is a piece-wise defined curve such that in each interval it is the linear combination of certain linearly independent functions $B_1, B_2, ... ,B_k$
- **Formal Definition**
    1. Let $J=[0, T]$ and consider the partition of  $J$ given by  $N + 1$ points $t_i\in J$, i.e. $I_1, I_2, ... ,I_N$ with $I_i=[t_i, t_{i + 1})$, $t_0=0$ and $t_{N + 1}=T$
    2. Let $I_0=[-1,1]$ and $B_1, B_2, ... ,B_k$ be $k$ linearly independent functions $B_i:I_0\longrightarrow \mathbb{R}$.
    3. Let $s_i:I_i\longrightarrow I_0$ given by

$$
    s_i(t)= 2\frac{t-t_i}{t_{i + 1}-t_i} - 1
$$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;iv. Let $\mathbf{y}_i^j \in\mathbb{R}^n$.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;v. A **generalized spline** from $J$ into $\mathbb{R}^n$ is a curve given by

$$
    f_r= (y_{i}^j)^\top \mathbf{B} \circ s_j(t) \text{ if } t\in I_j\ \ \ \ \ \ \ \ \ \ \ \ \ (\star)
$$

where $\mathbf{B}$ is the vector constructed by stacking the basis in a colunm vector

# Motivation

The term 'generalized splines' is borrowed from E. W. Cheney's book, *Introduction to Approximation Theory*, where he introduces generalized polynomials.
A generalized spline is to a generalized polynomial what a spline is to a polynomial.
Essentially, a gspline is a continuous piecewise function that passes through a sequence of points $\(\{q_1, \ldots, q_{N+1}\}\)$ at times $\(\{t_1, \ldots, t_{N+1}\}\)$.
Between each time interval $\([t_i, t_{i+1}]\)$, there is an arbitrary function $\(f_i: [t_i, t_{i+1}] \to \mathbb{R}^n\)$ that defines the gspline on $\([t_i, t_{i+1}]\)$.


Generalized splines appear naturally in trajectory optimization problems when waypoint constraints without prescribed times are added.
In other words, if we aim to optimize a motion that passes through a sequence of positions without specifying the times at which the motion meets each point, we will encounter generalized splines.

Generalized splines trajectories arising in such kind of optimization problems are uniquely characterized by the waypoints that they attain, the time at which each waypoint is meet and the speed, and possible higher order derivatives at the boundaries of the time interval, i.e, $t_0$ and $t_{N+1}$.
Moreover, such a relation is synthesised in the expression of the type

$$
\mathbf{A}(\boldsymbol\tau)\mathbf{y} = \mathbf{b}(\mathbf{w})\ \ \ \ \ \ \ \ \ \ \ \ \ \ (0)
$$

where $\boldsymbol{\tau} \in \mathbb{R}^N$ is a column vector representing the lenght of the time interval between consecutive waypoints, $\tau_i = t_{i+1} - t_i$. $\mathbf{A}(\boldsymbol{\tau})$ is a matrix dependent on the time intervals $\boldsymbol{\tau}$, $\mathbf{b}(\mathbf{w})$ is a column vector dependent on the waypoints $\mathbf{w}$, the speed, and possibly higher-order derivatives at the boundaries, and $\mathbf{y}$ is a column vector that uniquely represents the curve at each interval.

The main challenge to build into a computer a trajectory optimization problems with waypoint constraints is to compute the derivatives of $\mathbf{y}$ with respect to $\mathbf{\tau}$.

This library provides a uniform and simple interface to formulate gradient-based optimization problems for waypoint-constrained trajectory planning. The library leverage on the representation (0) to compute the "derivatives of the splines" with respect to the time intervals (and possibly the waypoints) as the corresponding derivatives of $\mathbf{y}$


# Background

This library is aimed to find a trajectory passing trough a sequence of waypoints $\mathbf{w}\_0, ...,\mathbf{w}_{N + 1}$ such that the following integral is minized

$$
\Large I=\int_0^T \alpha_1\left\|\frac{\mathsf{d}\mathbf{q}}{\mathsf{d} t }\right\|^2 + \alpha_2 \left\|\frac{\mathsf{d}^2\mathbf{q}}{\mathsf{d} t^2 }\right\|^2 + \alpha_3\left\|\frac{\mathsf{d}^3\mathbf{q}}{\mathsf{d} t^3 }\right\|^2 +  \alpha_4\left\|\frac{\mathsf{d}^4\mathbf{q}}{\mathsf{d} t^4 }\right\|^2 \mathsf{d} t \ \ \ \ \ (1)
$$

It may be proven that such a problem can be subdivided in two steps

 1. Find the family of optimal curves that joint waypoints
 2. Compute time instants $\{t_0, t_1,  ...,t_N, t_{N + 1}\}$ where the optimal curves must be attached

The step 1. is done by solving a linear ordinary differential equation, regarless the fact that the times $t_i$ are unknown. One method to achieve 2. is to formulate an optimization problem (e.g. a gradient based one) for which decission variable is the times at which each waypoint is acchived.

## Optimal curves
We underline that this library leverages on the [general theory of linear ODEs](https://en.wikipedia.org/wiki/Linear_differential_equation).
It may be proven that any optimal of (1) solves the following linear ODE, which turn out to be the Euler-Lagrange equations at each interval $[t_i, t_{i + 1}]$

$$
-\alpha_1\frac{\mathsf{d}^2\mathbf{q}}{\mathsf{d} t^2 } + \alpha_2 \frac{\mathsf{d}^4\mathbf{q}}{\mathsf{d} t^4 } - \alpha_3\frac{\mathsf{d}^6\mathbf{q}}{\mathsf{d} t^6 } +  \alpha_4 \frac{\mathsf{d}^8\mathbf{q}}{\mathsf{d} t^8 } = 0\ \ \ \ \ (2)
$$

with the following boundary conditions

$$
\mathbf{q}(t_i) = \mathbf{w}\_i\ \ \ \ \ \ \ \mathbf{q}(t_{i+1}) = \mathbf{w}_{i+1}\ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)
$$

Because the ODE (2) is linear, we can compute its general suction depending on the value of the coefficients $\alpha_i$.

In fact, the general solution of (2) may be written as a piecewise function defined at each interval as

$$
\mathbf{q} = \sum_{i=1}^{n_b} \mathbf{y}\_i^j B_i(t) \ \ \ \ \text{if}\ \ \ \ t \in [t_{j}, t_{j + 1}]\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (4)
$$

where $n_b$, $B_i(t)$ depend on the coefficients $\alpha_i$ and $\mathbf{y}\_i^j$ are column vectors in which represents the curve uniquely at the interval $[t_j, t_{j + 1}]$.

If we stack the column vectors $\mathbf{y}\_i^j$ in a suitable way we obtain the column vector $\mathbf{y}$ used in (0). In fact, the equation (0) is obtained after applying to (4) the waypoint constrains and the boundary conditions.

After substituting (4) in (1) we obtain the following expression

$$
I=\mathbf{y}^\top \mathbf{Q}(\boldsymbol{\tau}) \mathbf{y}\ \ \ \ \ \ \ \ \ \ \ \ \ (5)
$$

Finally we can substitute (0) in (5) to obtain the representation of (1) subject to the waypoint constrains as a function of $N$ real variables:

$$
I=I(\boldsymbol{\tau})=\mathbf{b}^{\top}\mathbf{A}^{-\top}(\boldsymbol{\tau})\mathbf{Q}(\boldsymbol{\tau}) \mathbf{A}^{-1}(\boldsymbol{\tau})\mathbf{b}\ \ \ \ \ \ \ \ \ \ \ \ \ (6)
$$

# Software architecture

From the formalization of the optimization problem, we derive that a flexible and uniform methodology for the construction of the problem of optimizing (1) consists in designing an abstract representation of the basis $B_i(t)$ in (4) capable of building in an automatic fashion the constraint (0), the cost function (5) and their derivatives.

In fact, note that the input of any gradient-based optimizer is the expression (6)  and its derivatives.
This library provides a template class to represent the basis $B_i(t)$ and a series of procedures which utilizes these basis as an input and then generate (0), (6) and their derivatives.

This library provides the class `gsplines::basis::Basis` which represent an arbitrary set of linearly independent functions and the class `gsplines::GSpline` which implements the definition $\star$.
In addition this library provides the class `gsplines::functions::FunctionBase` and `gsplines::functions::FunctionExpression` that allows to define arbitrary curves in $\mathbb{R}^n$ with the possibility of performing algebraic operations and derivation.

## Arbitrary function definition and algebraic operations

The abstract class [`gsplines::functions::FunctionBase`](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/master/include/gsplines/Functions/FunctionBase.hpp) represents objects that can be evaluated at points of some interval of $\mathbb{R}$. To inherit from this class the auxiliary template class [`gsplines::functions::FunctionInheritanceHelper`](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/master/include/gsplines/Functions/FunctionInheritanceHelper.hpp) allows to define a custom function which $k$-derivative is implemented as another class using the [Curiously recurring template pattern](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern). For example, the exponential and other elementary functions are [declared here](https://github.com/rafaelrojasmiliani/gsplines_cpp/blob/fa7395af1719a9c560eba15ff4bc68584da2c7e7/include/gsplines/Functions/ElementalFunctions.hpp#L82).


The class `gsplines::functions::FunctionExpression` inherits from `gsplines::functions::FunctionBase` and contain an array of pointers to `gsplines::functions::FunctionBase` called `function_array_`. The evaluation operator of `gsplines::functions::FunctionExpression` evaluate the desired algebraic operation between the functions in `function_array_`.
This architecture allow to implement complex algebraic expressions by recursively calling the evaluation of each function in `function_array_`.


# Requirements

- `numpy`
- `scipy`
- `matplotlib`
- `coinor-libipopt-dev`
- `ros-noetic-ifopt`
- `libeigen3-dev`

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

- Rojas, Rafael A., Andrea Giusti, and Renato Vidoni. "Online Computation of Time-Optimization-Based, Smooth and Path-Consistent Stop Trajectories for Robots." Robotics 11.4 (2022): 70.
```
@article{rojas2022online,
  title={Online Computation of Time-Optimization-Based, Smooth and Path-Consistent Stop Trajectories for Robots},
  author={Rojas, Rafael A and Giusti, Andrea and Vidoni, Renato},
  journal={Robotics},
  volume={11},
  number={4},
  pages={70},
  year={2022},
  publisher={MDPI}
}
```
