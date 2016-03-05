disp 'loading data...'
cd('../')
addpath('utils/')
load('results/count_bigrams.mat')
Ts = {};
for i = 1:length(bins)-1
   Ts{i} = freqs(2*i-1:2*i,:);
end
sample = 1;
depth = 10;
w = [0.15, 0.3, 0.55]';
dataset = sprintf('Noise_depth_%d_sample_%d.txt',depth,sample);
[x0, pos, ~, ~] = load_noninterpolated_Immuno(dataset);

d = size(w,1);
t = size(bins,2);
noise = 0.001;
disp 'running mex...'
n = length(x0);
tic
Z1 = decode_state(viterbiMex(x0,Ts,double(bins),pos,w,noise,depth,d,t,n),3);
toc
disp 'running matlab...'
tic
Z0 = decode_state(viterbi(x0,Ts,bins,pos,d,w,noise,depth),3);
toc


sum(sum(abs(Z0-Z1)))
cd('utils/')