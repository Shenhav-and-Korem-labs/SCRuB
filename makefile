all: _site.yml
	Rscript -e "rmarkdown::render('index.Rmd')"
	Rscript -e "rmarkdown::render('R-tutorial.Rmd')"
	Rscript -e "rmarkdown::render('prescrub_intro.Rmd')"
	Rscript -e "rmarkdown::render('QIIME2-Tutorial.Rmd')"
	Rscript -e "rmarkdown::render('qiime2-intro.Rmd')"
	Rscript -e "rmarkdown::render('R-inputs.Rmd')"
