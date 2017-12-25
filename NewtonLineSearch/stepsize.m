function [alpha] = stepsize(f, g, x,rho, p,alpha_init, c )
% Calculates the alpha, or stepsize for the Newton method with backtracking
% line search.
%
% f - our objective function
% g - the calculated gradient vector at our point
% x - initial x vector
% p - search direction
% rho - constant between [0,1]. The proportion which we want to decrease our
% step size by.
% c - constant between [0,1]. Usually small, 10e-3 or smaller.
% alpha_init - initial step length
% alpha - returns the calculated step size.
alpha = alpha_init;
prev_f = f(x);
prev_x = x;
x =  prev_x + alpha*p;
curr_f = f(x);

while curr_f > prev_f + c*alpha*(g'*p)
    alpha = alpha*rho;
    
    x = prev_x + alpha*p;
    curr_f = f(x);
end

end

