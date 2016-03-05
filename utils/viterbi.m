function [S_star] = viterbi(x0,Ts,bins,pos,d,w,noise,depth)
D = 2^d;
n = length(x0);
V = zeros(n,D) - Inf;
% Fill in first row
for i = 1:D,
    state = de2bi(i-1,d);
    if noise == 0,
        V(1,i) = d*log(1/2) + logbinomialpdf(depth*x0(1),depth,state*w);
    else
        sigma = max(sqrt(depth*state*w*(1-state*w)*noise),1e-3);
        V(1,i) = d*log(1/2) + log(normpdf(depth*x0(1),depth*state*w,sigma));
    end
end
for i = 2:n,
    dist = pos(i)-pos(i-1);
    if dist <= 0 % gene restarts
        T = ones(2,2)/2;
    else
        [~,idx] = max(bins>=dist);
        T = Ts{idx-1};
    end
    for j = 1:D,
        for prev = 1:D,
            state = de2bi(j-1,d);
            bin_prev = de2bi(prev-1,d);

            % Probability of emission given state = bin
            if noise == 0,
                log_em_prob = logbinomialpdf(depth*x0(i),depth,state*w);
            else
                sigma = max(sqrt(depth*state*w*(1-state*w)*noise),1e-3);
                log_em_prob = log(normpdf(depth*x0(i),depth*state*w,sigma));
            end

            % Transition probability from previous state to here
            log_tr_prob = log(1);
            for k = 1:d,
                log_tr_prob = log_tr_prob + log(T(bin_prev(k)+1,state(k)+1));
            end

            V(i,j) = max(V(i,j),log_em_prob + log_tr_prob + V(i-1,prev));
        end
    end
end
for j = 1:n,
	[~,i] = max(V(j,:));
	S_star(j) = i;
end
end
