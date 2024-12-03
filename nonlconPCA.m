function [c,ceq] = nonlconPCA(x, tgt, fl, cv, D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = ((fl' * x)' * cv * (fl'*x) + x' * D * x) - tgt^2;
ceq = [];
end

