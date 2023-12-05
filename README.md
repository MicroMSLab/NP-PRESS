# NP-PRESS

This repository consists of two parts, FUNEL and simRank. We offer an online service for FUNEL and simRank, accessible at [https://npcompass.zju.edu.cn/](https://npcompass.zju.edu.cn/). If you prefer to use the codes locally, please follow the instructions below.

## FUNEL

FUNEL is written in R and is used to identify and annotate features from $MS^2$ data, and search compounds in a database. or more details, you can refer to the [README file](./FUNEL/README.md) in the project.

## simRank

SimRank is written in Python and has two modules:

1.	simRank-Filtering: filter features in samples based on their $MS^2$ spetra simRank similarity with the spectra in reference/control samples.
2.	simRank-Networking: cluster $MS^2$ spectra in samples based on simRank similarity scores.


For more details, you can refer to the [README file](./simRank/README.md) in the project.

## About Environment

While both FUNEL and simRank make use of the conda environment, we highly recommend employing distinct environments for each. Integrating Python and R within the same environment may not be the most optimal choice, and we have not tested these codes in such a configuration. However, if you are interested or have a need to do so, we welcome Pull requests.

Attention: Using Miniforge or other open-source conda distributions may lead to unexpected errors. We highly recommend using Anaconda for a more reliable experience.

## Citation

If you use this tool, please cite it by referencing the following publication:
[Author(s). Title of the Work. Journal/Conference, Year. DOI or URL]