function [x0, pos, gmints, gmnames] = load_noninterpolated_Immuno(dataset, no_genes=50)

raw = load(dataset);

gmints = raw(:,1);
pos = raw(:, 2) + 5000;
x0 = raw(:, 3);

% only take first `no_genes` genes
if nargin > 1
  % don't know where the gene numbers are starting, so shift them to 1.
  mask = (gmints - min(gmints) + 1) < no_genes;
  x0 = x0(mask);
  pos = pos(mask);
end
end

##
## example:
##
## load_noninterpolated_Immuno("path/to/file.cleaned.txt");
##
