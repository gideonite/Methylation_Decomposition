function [Z0] = load_truth_Immuno(dataset, no_genes)

raw = load(dataset);
gmints = raw(:,1);
pos = raw(:,2);
Z0 = raw(:,3:5);

if nargin > 1
  % don't know where the gene numbers are starting, so shift them to 1.
  mask = (gmints - min(gmints) + 1) < no_genes;
  Z0 = Z0(mask,:);
  pos = pos(mask,:);
end
