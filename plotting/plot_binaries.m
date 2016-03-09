addpath('../utils/');

## filenames are arrays of the form ["path" "to" "file" "filename"]
function plot_binary(noninterpolated_filename, results_filename,
                     truth_filename, no_genes=50) 
  load(results_filename);
  [x0, pos, ~, ~] = load_noninterpolated_Immuno(noninterpolated_filename, no_genes);

  pos = pos - 5000;
  ## double check I've sorted
  [w, I] = sort(w);
  Z = Z(:, I);
  keyboard
  Z0= Z0(1:size(Z,1),1:size(Z,2));
  assert(isequal(size(Z),size(Z0)))
  diff_ratio = sum(sum(abs(Z-Z0))) / prod(size(Z))

  
end

##
## example
##

plot_binary('/home/gideon/Data/methylation-decomposition/ImmunoMix3/Noise_depth_20_sample_8.cleaned.txt', '/home/gideon/Data/methylation-decomposition-results/hmm3.0_Immuno/20_8.mat') 

## pos = pos - 5000;
## % double check I'v sorted
## [w, I] = sort(w);
## Z = Z(:, I);
## Z0= Z0(1:size(Z,1),1:size(Z,2));
## assert(isequal(size(Z),size(Z0)))
## diff_ratio = sum(sum(abs(Z-Z0))) / prod(size(Z))
## % Plot the true vs. estimate binaries of the first gene
## d = length(w);
## clf;
## gmints(length(gmints)+1)=length(x0)+1;
## for c = 1:3
##     figure;
##     for i = 1:d
##         subplot(d,1,i);
##         plot(pos(gmints(c):gmints(c+1)-1), Z0(gmints(c):gmints(c+1)-1, i),'ro');
##         hold on
##         plot(pos(gmints(c):gmints(c+1)-1), Z(gmints(c):gmints(c+1)-1, i),'b*');
## %         legend('truth','estimated', 'Location','southeast');
## %         title(sprintf('binary curve %d with weight %0.2f', i, num2str(w(i))));
##         axis([-5000,5000,-0.1,1.1]);
##     end
## end
