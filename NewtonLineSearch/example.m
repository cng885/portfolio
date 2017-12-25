
% Example File
%
% @author Chase Ginther
% @date 2017.12.11

syms x

%%Our objective function
f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

%%Inline function represenitng our gradient. 
gradient = @(x) [400*x(1)^3-2*x(1)*(200*x(2)-1)-2 ;
                    200*x(2)-200*x(1)^2];

%%Inline function to calculate the hessian
hessian = @(x) [1200*x(1)^2-2*(200*x(2)-1) -400*x(1); 
                    -400*x(1) 200];
           



c = 1e-3;
rho = 0.2;
pointa = [1.2; 1];
pointb = [-1.2;1];

[min_a, table_a] = newtonlinesearch(1, f, pointa, hessian, gradient, c, rho);
table_a
[min_b, table_b] = newtonlinesearch(1, f, pointb, hessian, gradient, c, rho);
table_b