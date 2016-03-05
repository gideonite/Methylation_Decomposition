function [x0, pos, gmints, gmnames] = load_noninterpolated_Immuno(dataset, k)
% One of those two delimiters should work
delimiterIn = '\t';
% delimiterIn = ' ';

datapath = '/research-projects/mlgroup2/yusun/MCD/data/ImmunoMix3/';
headerlinesIn = 0;
raw = importdata([datapath dataset],delimiterIn, headerlinesIn);
x0 = raw.data(:, 2);
pos = raw.data(:, 1) + 5000;
[gmnames, gmints]= unique(raw.textdata,'stable');
if nargin > 1
	% only take first k genes
	x0 = x0(1:gmints(k+1)-1);
	pos = pos(1:gmints(k+1)-1);
	gmnames = gmnames(1:k);
	gmints = gmints(1:k);
end
end
