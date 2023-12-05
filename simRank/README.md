# simRank

This part contains the Python code for simRank. SimRank is an algorithm used to measure the similarity between two given $MS^2$ spectra. It includes two modules: `simRank-Filtering` and `simRank-Networking`, which will be described in the following sections.

## Init

We recommend using conda to manage the version of `Python` and its packages. For your convenience, we provide the `env.yaml` file, which can be used to create the required environment.

```bash
conda env create -f env.yaml
```

## simRank-Filtering

This script is designed to filter features in the experimental sample. It does this by comparing the $MS^2$ spectra of the sample with those in reference/control samples, using simRank similarity as the basis for comparison. It outputs a list of features after filtering. 

The following parameters are required:

`sample_path` is the path for experimental data file (.mzXML)

`peaktable_path` is the path for feature list before filtering, can be the output from FUNEL

`control_dir_list` is the list of folders containing reference/control data files

`output_dir` is the path for the output file

`scorematrix_path` is the path for the probability matrix

`inten_thresh` is the intensity cutoff for peaks in normalized $MS^2$ spectra

`rt` is the retention time tolerance to merge spectra from a same percursor

`ppm` is the m/z tolerance

`pm_tolerance` is the precursor difference allowed to compare two spectra

`msdelta` is the m/z tolerance for two precursors to be considered as the same

`ms2delta` is the m/z tolerance for two fragments in an $MS^2$ spectrum to be considered as the same

`min_mz_num` is the minimum number of fragments required in a valid $MS^2$ spectrum

`if_merge_database`, if True, spectra from the same precursor in each reference/control sample will be merged

`if_merge_database_byenergy`, if False, spectra from the same precursor in each reference/control sample will all be merged even with different collision energies

`if_merge_samples_byenergy`, if False, spectra from the same precursor in the experimental sample will all be merged even with different collision energies

`remove_precursor`, if True, fragments within 17 Da of the precursor ion will be removed in an $MS^2$ spectrum

`simrank_thr` is the simRank score cutoff for feature filtering (recommend to use a value above 16)


## simRank-Networking

This script is designed to perform clustering analysis of $MS^2$ spectra from the experimental sample and reference/control samples. It outputs a GraphML file that can be visualized by Cytoscape.The following parameters are the same as those for simRank-Filtering. 