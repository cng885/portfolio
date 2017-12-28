function [ w ] = rungeKutta( f, alpha, h, t, N)
%This function is for the Runge Katta method.
% This is a numerical approximation solver for an ODE
%
%@param f the function of the ODE
%@param alpha our initial condition
%@param h specified step size
%@param t time array
%@param N array size
%
%@author Chase Ginther
%@date 2016.11.13


w = zeros(size(t)-1);
w(1) = alpha;

for i = 1:N
    k1 = h*(f(t(i), w(i)));
    k2 = h*(f(t(i)+h/2, w(i)+1/2*k1));
    k3 = h*(f(t(i)+h/2, w(i)+1/2*k2));
    k4 = h*(f(t(i)+1, w(i)+k3));

    w(i+1) = w(i)+1/6*(k1+2*k2+2*k3+k4);
end

end
