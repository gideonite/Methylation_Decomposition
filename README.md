# Methylation_Decomposition
Decomposing methylation curves from a mixture of cell population

## Getting the Dataset to Work
This demo uses a dataset called ImmunoMix3, which contains many simulated binary methylation curves for three cell populations as ground truth, with the observed curve at various levels of noise (i.e. depth of coverage). The dataset is available at <br/>
https://www.dropbox.com/s/ovksanvsgrmmg17/ImmunoMix3.zip?dl=0 <br/>
After you download the dataset, please set the datapath variable at those two file loaders: <br/>
- <code> utils/load_truth_Immuno.m </code> <br/>
- <code> utils/load_noninterpolated_Immuno.m </code>


## Compiling the mex function
By default, test_HMM.m uses the very fast mex function that we wrote for viterbi, which performs the same functionality as viterbi.m but is about 1000x faster (as timed by vitFullTest.m). Starting to use the mex is very simple, just do: <br/>
<code> mex utils/viterbiMex.cpp </code> <br/>
To use the slow MATLAB function instead, please do/undo comment for the corresponding lines in test_HMM.m.