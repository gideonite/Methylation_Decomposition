function [ w, fval ] = quadmin(Z, x, opts)
d = size(Z, 2);
[w, fval] = quadprog(double(Z)'*Z, -Z'*x, [],[], ...
    ones(1,d),1, ...
    zeros(d,1),[], [], opts);
end

% quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opts) 