function [ output ] = logbinomialpdf(k,n,p)
output = sum(log(n-k+1:n)) - sum(log(1:k))...
    + k*log(p) + (n-k)*log(1-p);
end

