function [Z,w,curve_errors] = test_HMM(opts,Ts,bins,t,x0,pos,iters,noises,depth)
n = size(x0, 1);
w = rand(t, 1); w = w./sum(w);
curve_errors = zeros(sum(iters),1);
for r = 1 : length(iters)
    for i = 1 : iters(r)
        idx = sum(iters(1:r-1))+i;
        % Z = decode_state(viterbi(x0,Ts,bins,pos,t,w,noises(r),depth),t);			% this does NOT use the mex function
        Z = decode_state(viterbiMex(x0,Ts,double(bins),pos,w,noises(r),depth,t,size(bins,2),length(x0)),t);
        w = quadmin(Z, x0, opts);
        curve_errors(idx) = sqrt(sum((x0 - Z*w).^2)/n);
    end
end
[w, I] = sort(w);
Z = Z(:, I);
end
