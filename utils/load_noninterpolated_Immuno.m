function [x0, pos, gmints, gmnames] = load_noninterpolated_Immuno(dataset, no_genes)

raw = load(dataset);

geneno = raw(:,1);
pos = raw(:, 2) + 5000;
x0 = raw(:, 3);

[gmnames, gmints] = unique(geneno);

% only take first `no_genes` genes
if nargin > 1
  % don't know where the gene numbers are starting, so shift them to 1.
  mask = (geneno - min(geneno) + 1) < no_genes;

  x0 = x0(mask);
  pos = pos(mask);

  gmnames = gmnames(1:no_genes);
  gmints = gmints(1:no_genes);
end

end


%%
%% example:
%%
%% load_noninterpolated_Immuno("Noise_depth_20_sample_8.cleaned.txt");
%%
