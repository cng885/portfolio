function [ w ] = euler( f, h, i, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



t = I(1):h:I(2);

w = zeros(size(t)-1);
w(1) = alpha;

for i = 1:N-1
    w(i+1) = w(i)+h*f(t(i),w(i));
end

end

