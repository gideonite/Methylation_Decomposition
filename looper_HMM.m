function [] = looper_HMM()
Zerrs = ones(10,20);
werrs = ones(10,20);
w0 = [0.15,0.3,0.55]';
addpath('utils/')
load('results/count_bigrams.mat')
Ts = {};
for i = 1:length(bins)-1
    Ts{i} = freqs(2*i-1:2*i,:);
end
opts = optimoptions('quadprog','Display','off',...
    'Algorithm', 'interior-point-convex' );
k = 50;
t = 3;
iters = [2,2,2];
noises = [2,1.5,1.0];

datadir = '/scratch/gmd87/ImmunoMix3/';

for sample = 1:20
    disp(sample);
    truth_dataset = [datadir sprintf('Solutions_OneThousandGene_%d.cleaned.txt',sample)];
    Z0 = load_truth_Immuno(truth_dataset, k);
    for depth = (1:10) * 10
        dataset = [datadir sprintf('Noise_depth_%d_sample_%d.cleaned.txt',depth,sample)];
        [x0, pos, ~, ~] = load_noninterpolated_Immuno(dataset,k);
        [Z,w,curve_errors] = test_HMM(opts,Ts,bins,t,x0,pos,iters,noises,depth);
        save(sprintf('results/hmm3.0_Immuno/%d_%d.mat',depth,sample),...
            'Z','w','curve_errors');
        %%%%%%
        Zerrs(depth/10,sample) = sum(sum(abs(Z-Z0))) / prod(size(Z));
        werrs(depth/10,sample) = norm(w-w0,2);
    end
end
save('results/hmm3.0_Immuno/errors.mat', 'Zerrs', 'werrs');

end
