for rmd in  annotation-download-format.Rmd wgs-lowmap-bias.Rmd wgs-coverage-tracks.Rmd PopSV-methodBenchmark-twin-reliability.Rmd PopSV-methodBenchmark-twin.Rmd PopSV-methodBenchmark-cagekid.Rmd PopSV-catalog-overview.Rmd PopSV-NA12878-assembly-lowMap.Rmd PopSV-repeatEnr.Rmd PopSV-functionalImpact.Rmd
do
    Rscript -e "rmarkdown::render('src/$rmd')"
    mv -f src/${rmd%Rmd}md reports/
    rm -rf reports/${rmd%.Rmd}_files
    mv src/${rmd%.Rmd}_files reports/
done
