This repository hosts the source code and instructions to reproduce the analysis and results from our study on low-mappability regions and copy number variation.

If you just want the CNV catalog, find them in the [FigShare repository](https://figshare.com/s/8fd3007ebb0fbad09b6d), or directly [here](https://ndownloader.figshare.com/files/3638574).

**To review the code and resulting graphs/numbers, have a look at the R-markdown reports** in the [`reports`](https://github.com/jmonlong/reppopsv/tree/master/reports) folder and scripts in the  [`src`](https://github.com/jmonlong/reppopsv/tree/master/src) folder.

To rerun the analysis, follow these steps:

# Step 1: Download the relevant data

The necessary data has been deposited on [FigShare](https://figshare.com/s/8fd3007ebb0fbad09b6d). Depending on the analysis, you might not need to download all the data.
Still, the easiest way is to download all the data and unzip it in the `data` folder.

```sh
mkdir -p data
cd data
wget https://ndownloader.figshare.com/articles/2007630?private_link=8fd3007ebb0fbad09b6d -o figshare.zip
unzip figshare
tar -xzvf PopSV-NA12878-lowmap.tar.gz
```

# Step 2: Install R dependencies

Many different packages are used throughout the analysis. The commands to install them are written in the `installDependencies.R`. 
To install all the necessary packages open R and run `source("installDependencies.R")`.

# Step 3: Compile the R-markdown reports

The raw R-markdown reports are located in the `src` folder. To recompile them simply run:

```r
library(rmarkdown)
render("XXX.Rmd")
```

Of note, the `annotation-download-format.Rmd` report should be compiled before anything because it downloads and format public annotations.

Or, you can compile all the reports using:

```sh
sh compileAll.sh
```

You can already see the reports produced by these scripts in the `reports` folder.

# Notes

The code was tested on fresh dockerized Ubuntu with [R 3.3.1](https://github.com/jmonlong/myDockerfiles/blob/reppopsv/R.3.3.1/Dockerfile). Windows is not recommended as it doesn't support the `parallel` package.
