# Newton Line Search Routine

A local optimization method for finding a local minima. Adopted from Numerical Optimization by Nocedal and Wright. This contains an implementation in Matlab. Our objective function must be twice continuously differentiable. 

## Files

### NewtonLineSearch.m
Our newton line search routine.  The first step of this routine is that our Hessian must be symmetric positive definite. We can use a modified Hessian method to ensure this iwthout disturbing our solution. We simply add a multiple of the smallest eigenvalue of the Hessian to ensure that it is positive. We then must solve for our search direction, denoted $p$ by solving $Bp=-g$. We then use this search direction and call our step size routine to generate our step size we will take. 

### stepsize.m
Calculates the step size we will take. Iterates until we have researched the optimum step size to take.

### example.m

Example problem when we find the minima of a 2 dimensional Rosenbrock function. 
