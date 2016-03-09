addpath('utils/')

filename = 'jerry2';
load(['data/noninterpolated/' filename '_truth.mat']);
load(['results/enter file name here' filename '.mat']);
[x0, pos, gmints, gmnames] = load_noninterpolated([filename '.dat'],3);
pos = pos - 5000;
% double check I'v sorted
[w, I] = sort(w);
Z = Z(:, I);
Z0= Z0(1:size(Z,1),1:size(Z,2));
assert(isequal(size(Z),size(Z0)))
diff_ratio = sum(sum(abs(Z-Z0))) / prod(size(Z))
% Plot the true vs. estimate binaries of the first gene
d = length(w);
clf;
gmints(length(gmints)+1)=length(x0)+1;
for c = 1:3
    figure;
    for i = 1:d
        subplot(d,1,i);
        plot(pos(gmints(c):gmints(c+1)-1), Z0(gmints(c):gmints(c+1)-1, i),'ro');
        hold on
        plot(pos(gmints(c):gmints(c+1)-1), Z(gmints(c):gmints(c+1)-1, i),'b*');
%         legend('truth','estimated', 'Location','southeast');
%         title(sprintf('binary curve %d with weight %0.2f', i, num2str(w(i))));
        axis([-5000,5000,-0.1,1.1]);
    end
end
