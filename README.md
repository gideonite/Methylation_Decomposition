# Methylation_Decomposition
Decomposing methylation curves from a mixture of cell population

## Getting the Dataset to Work
This demo uses a dataset called ImmunoMix3, which contains many simulated binary methylation curves for three cell populations as ground truth, with the observed curve at various levels of noise (depth of coverage). The dataset is available at <br/>
https://www.dropbox.com/sh/8jo1g4fn67752f9/AAAqk11QLx49Ld5gGt0AcQaQa?dl=0 <br/>
After you download the dataset, please set the datapath variable at the two file loaders: <br/>
utils/load_truth_Immuno.m <br/>
utils/load_noninterpolated_Immuno.m <br/>


## Compiling the mex function
By default, test_HMM.m uses the very fast mex function that we wrote for viterbi, which performs the same functionality as viterbi.m but is about 1000x faster (as timed by vitFullTest.m). Starting to use the mex is very simple, just do: <br/>
mex utils/viterbiMex.cpp <br/>
To use the slow MATLAB function instead, please do/undo comment for the corresponding lines in test_HMM.m.