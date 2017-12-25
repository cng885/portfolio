function [ x_new, table_result ] = newtonlinesearch( alpha_init, f,x_0, hessian, gradient, c, rho )
% Implementation of a Newton (back-tracking) Line Search. Uses a seperate
% function to calculate the step size.
%
% alpha_init - initial step length
% f - the function handle
% x_0 - initial x vector of dimension n.
% gradient - the gradient function handle
% hessian - the hessian function handle 
% rho - constant between [0,1]. The amount we will decrease each step size
% by.
% c - constant between [0,1]. Usually smaller than 10e-3.
% 
% @author Chase Ginther
% @date 2017.12.05

max_iter = 100;
x_new = x_0 ; 
iter = zeros(max_iter,1);
alphas = zeros(max_iter, 1);
gradient_norm = zeros(max_iter, 1);
fval = zeros(max_iter, 1);

for i=1:max_iter
    
    x_old = x_new;
    h = hessian(x_old);
    g = gradient(x_old);
    
    %get the smallest eigenvalue from the hessian for Hessian Modification
    lambdas = eig(h);
    if (min(lambdas) < 0)
        Bk = h - 2*min(lambdas)*eye(size(h));
    else
        Bk = h;     
    end 
 
    %solve the negative gradient for the search direction
    pk = (-Bk\g);
    
    
    %use the backtrack line search to calc the step size
    alpha = linesearch(f, g, x_old, rho, pk, alpha_init, c);

    x_new = x_old +  alpha* pk;
    
    %populate our arrays to be stored in nice table. 
    iter(i) = i;
    alphas(i) = alpha;
    gradient_norm(i) = norm(g);
    fval(i) = f(x_new);
end


table_result = table(iter, alphas, gradient_norm, fval);

end

