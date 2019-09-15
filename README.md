# sib

System Identification Toolbox

## Description

This is a Matlab toolbox for parameter identification of dynamic systems.
The focus of the project is to obtain good models to the system, not fast algorithms.
The toolbox is also known as Slow is Better (sib).

Implemented methods:

* Stieglitz-McBride
* Predicition Error Method with ARX structure
* Predicition Error Method with ARMAX structure
* Predicition Error Method with ARMAX structure with improved convergence algorithm
* Predicition Error Method with OE structure
* Predicition Error Method with OE structure with improved convergence algorithm
* Predicition Error Method with BJ structure
* Predicition Error Method with BJ structure with improved convergence algorithm


## Install

Two steps are necessary to install **sib**:

1. Some algorithms are written in *C* and need to be compliled with mex.
Run the matlab file ./src/C/install.m.

2. Add the folder ./src/ to the *path* of Matlab.

## Use

Please check the *examples* folder.

#### Example of ARX model

```matlab
N = 1000;
u = square((1:N)*(2*pi)/100)';
G = tf([0 1],[1 -0.9],1);
H = tf([1 0],[1 -0.9],1);

y = lsim(G,u) + lsim(H,randn(N,1));

[theta, m] = sib_arx(u, y, 1, 1, 1)
```


## Contributors

Diego Eckhard - diegoeck@ufrgs.br - @diegoeck
