function [] = plot_curve_non()
filename = 'jerry2';
load(['results/hmm2.9_' filename '.mat']);
[x0, pos, gmints, gmnames] = load_noninterpolated([filename '.dat']);
plot_curve_non_inside(Z, w, x0, gmints, gmnames, pos)
end


function [] = plot_curve_non_inside(Z, w, x0, gmints, gmnames, pos)
% Plots the first n curves, estimated and truth
xe = Z*w;
pos = pos - 5000;
n = 4;
% Only if n is equal to the total # of genes
% gmints(n+1) = length(pos)+1;
clf;
fig = figure;
for i = 1:n
    subplot(n,1,i);
    plot(pos(gmints(i):gmints(i+1)-1), x0(gmints(i):gmints(i+1)-1),'r-o');
    hold on
    plot(pos(gmints(i):gmints(i+1)-1), xe(gmints(i):gmints(i+1)-1),'b-*'); 
    legend('truth','estimate', 'Location','southeast');
    title(gmnames(i));
    axis([-5000,5000,-0.1,1.1]);
end
end