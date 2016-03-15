% TODO: make this into some sort of proper function.

% TODO: create a file called "global_variables.m" and all this stuff
% there.
addpath('utils/');
datadir = '/scratch/gmd87/ImmunoMix3/';

no_genes = 50;

%% load the transition matrix from bigrams in the training data. See
%% `count_bigrams.py`.
load('results/count_bigrams.mat')
Ts = {};
for i = 1:length(bins)-1
    Ts{i} = freqs(2*i-1:2*i,:);
end
%%

%% set optimization options
opts = optimoptions('quadprog','Display','off',...
    'Algorithm', 'interior-point-convex' );
%% 

%% randomly initialize w but normalize s.t. sum(w) == 1.
t = 3; % TODO the number of binary methylation signal to learn.
w = rand(t, 1); w = w./sum(w);
%%

ground_truth_w = [0.15,0.3,0.55]';

samples = 1:20; %samples = 1:1;
sequencing_depths = 10:10:100; %sequencing_depths = 100;
noise_levels = 1:.1:2;
w_errors = zeros(length(samples), length(sequencing_depths), length(noise_levels));
for i = 1:length(samples)
  disp(['sample ' int2str(samples(i))]);
  truth_dataset = [datadir sprintf('Solutions_OneThousandGene_%d.cleaned.txt',samples(i))];

  Z0 = load_truth_Immuno(truth_dataset, no_genes);

  for j = 1:length(sequencing_depths)
    dataset = [datadir sprintf('Noise_depth_%d_sample_%d.cleaned.txt', ...
                               sequencing_depths(j), samples(i))];
    [x0, pos, ~, ~] = load_noninterpolated_Immuno(dataset, no_genes);

    iters = repmat(2,1,length(noise_levels)); % [2, 2, ..., 2]
    for k = 1 : length(noise_levels)
      avg_w = 0;
      for unused = 1 : iters(k)
        viterbiOutput = viterbiMex(x0, Ts, double(bins), pos, w, ...
                                   noise_levels(k), sequencing_depths(j), t, size(bins,2), length(x0));
        Z = decode_state(viterbiOutput,t);
        w = quadmin(Z, x0, opts);
        avg_w = avg_w + w;
      end

      avg_w =  avg_w ./ iters(k);
      w_err = norm(avg_w - ground_truth_w,2);
      w_errors(i, j, k) = w_err;
    end
  end
end

save('results/hmm3.0_Immuno/w_errors.mat', 'w_errors');

%% save('results/hmm3.0_Immuno/errors.mat', 'Zerrs', 'werrs');
