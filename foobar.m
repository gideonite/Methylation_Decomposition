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

%% samples = 1:20;
samples = 1:1;
for sample = samples
  disp('sample');
  disp(sample);
  truth_dataset = [datadir sprintf('Solutions_OneThousandGene_%d.cleaned.txt',sample)];

  Z0 = load_truth_Immuno(truth_dataset, no_genes);

  for sequencing_depth = 100
  %% for sequencing_depth = 10:10:100
    dataset = [datadir sprintf('Noise_depth_%d_sample_%d.cleaned.txt',sequencing_depth,sample)];
    [x0, pos, ~, ~] = load_noninterpolated_Immuno(dataset, no_genes);

    noises = 1:.1:2;
    iters = repmat(2,1,length(noises)); % [2, 2, ..., 2]
    for noise_level = 1 : length(noises)
      avg_w = 0;
      for unused = 1 : iters(noise_level)
        viterbiOutput = viterbiMex(x0, Ts, double(bins), pos, w, ...
                                   noises(noise_level), sequencing_depth, t, size(bins,2), length(x0));
        Z = decode_state(viterbiOutput,t);
        w = quadmin(Z, x0, opts);
        avg_w = avg_w + w;
      end

      avg_w =  avg_w ./ iters(noise_level);
      w_err = norm(avg_w - ground_truth_w,2);
    end
  end
end
%% save('results/hmm3.0_Immuno/errors.mat', 'Zerrs', 'werrs');
