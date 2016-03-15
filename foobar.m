% TODO: create a file called "global_variables.m" and dump all this
% sh**t there.
addpath('utils/');
datadir = '/scratch/gmd87/ImmunoMix3/';

k = 50;

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

%% samples = 1:20;
samples = 1:1;
for sample = samples
  disp('sample');
  disp(sample);
  truth_dataset = [datadir sprintf('Solutions_OneThousandGene_%d.cleaned.txt',sample)];

  Z0 = load_truth_Immuno(truth_dataset, k);

  for sequencing_depth = 100
  %% for sequencing_depth = 10:10:100
    dataset = [datadir sprintf('Noise_depth_%d_sample_%d.cleaned.txt',sequencing_depth,sample)];
    [x0, pos, ~, ~] = load_noninterpolated_Immuno(dataset,k);

    noises = 1:.1:2;
    iters = repmat(2,1,length(noises)); % [2, 2, ..., 2]
    for noise_level = 1 : length(noises)
      for unused = 1 : iters(r)
        viterbiOutput = viterbiMex(x0, Ts, double(bins), pos, w, ...
                                   noises(noise_level), sequencing_depth, t, size(bins,2), length(x0));
        Z = decode_state(viterbiOutput,t);
        w = quadmin(Z, x0, opts);
      end
    end

    %% [Z,w,curve_errors] = test_HMM(opts,Ts,bins,t,x0,pos,iters,noises,sequencing_depth);
    %% save(sprintf('results/hmm3.0_Immuno/%d_%d.mat',sequencing_depth,sample),...
    %%     'Z','w','curve_errors');
    %% %%%%%%
    %% Zerrs(sequencing_depth/10,sample) = sum(sum(abs(Z-Z0))) / prod(size(Z));
    %% werrs(sequencing_depth/10,sample) = norm(w-w0,2);
  end
end
%% save('results/hmm3.0_Immuno/errors.mat', 'Zerrs', 'werrs');
