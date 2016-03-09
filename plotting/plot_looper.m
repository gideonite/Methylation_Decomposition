addpath('utils/')
load('results/hmm3.0_Immuno/errors.mat', 'Zerrs', 'werrs');

figure
avgs = mean(Zerrs,2);
devs = std(Zerrs,0,2);
H = shadedErrorBar((1:10)*10,avgs,devs, '*b-');
ylabel('Fraction of error on Z')
xlabel('Depth of Coverage')

figure
avgs = mean(werrs,2);
devs = std(werrs,0,2);
G = shadedErrorBar((1:10)*10,avgs,devs, '*r-');
ylabel('Euclidean distance on w')
xlabel('Depth of Coverage')

% figure
% for i = 1:20
%     hold on
%     scatter((1:10)*10, Zerrs(:,i), 20, 'filled')
%     hold on
% end
% ylabel('Fraction of error on Z')
% xlabel('Depth of Coverage')
% 
% figure
% for i = 1:20
%     hold on
%     scatter((1:10)*10, werrs(:,i), 20, 'filled')
%     hold on
% end
% ylabel('Euclidean distance on w')
% xlabel('Depth of Coverage')