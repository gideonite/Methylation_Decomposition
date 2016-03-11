# Methylation_Decomposition
Decomposing methylation curves from a mixture of cell population

## Getting the Dataset to Work
This demo uses a dataset called ImmunoMix3, which contains many simulated binary methylation curves for three cell populations as ground truth, with the observed curve at various levels of noise (i.e. depth of coverage). The dataset is available at <br/>
https://www.dropbox.com/s/ovksanvsgrmmg17/ImmunoMix3.zip?dl=0 <br/>

After you download the dataset, symlink the data directory to `./data`. That is, `cd` into the repo and run

```
$ ln -s /path/to/data/ImmunoMix3 data
```

## Compiling the mex function

### Matlab
By default, test_HMM.m uses the very fast mex function that we wrote for viterbi, which performs the same functionality as viterbi.m but is about 1000x faster (as timed by vitFullTest.m). Starting to use the mex is very simple, just do: <br/>
<code> mex utils/viterbiMex.cpp </code> <br/>
To use the slow MATLAB function instead, please do/undo comment for the corresponding lines in test_HMM.m.

### Octave

For more info, see [this](https://www.gnu.org/software/octave/doc/interpreter/Getting-Started-with-Mex_002dFiles.html#Getting-Started-with-Mex_002dFiles)

#### tl;dr

```sh
mkoctfile --mex utils/viterbiMex.cpp
```

## Run it

```sh
$ matlab
>> looper_HMM
     1
     2
     ...
```
You should see some numbers coming out as the looper does its thing. The results are stored in `results/hmm3.0_Immuno/<depth>_<sample>`.

# TODO

- document the input data and the output data as the corresponding data cleaning script(s) for converting various formats to the standard input.
