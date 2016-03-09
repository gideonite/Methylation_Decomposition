for depth = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
  for sample = 1:10
    load(strcat("~/Data/methylation-decomposition-results/hmm3.0_Immuno/",
                int2str(depth), "_", int2str(sample), ".mat"));

    for zcol = 1:3
      subplot(5,1,zcol+1);
      plot(Z(:,zcol),'-');
      ylim([-0.05,1.05]);
      title(strcat("decomposed signal no. ", int2str(zcol)));
    end

    X = dlmread(strcat('/home/gideon/Data/methylation-decomposition/ImmunoMix3/Noise_depth_',int2str(depth), '_sample_', int2str(sample), '.txt'), '\t');
    wndw = size(X,1) / size(Z,1);

    mask = zeros([size(X,1),1]);
    step = floor(size(X,1) / size(Z,1));
    for i = 1:size(X,1)
      if mod(i,step) == 0
        mask(i) = 1;
      end
    end
    mask = logical(mask);

    X = X(:,3);
    X = X(mask);

    subplot(5,1,5);
    plot(X,'-');
    ylim([-0.05,1.05]);
    title("original signal");

    print("pdf", strcat("plots/", int2str(depth), "_", int2str(sample), ".pdf"));

  end
end


