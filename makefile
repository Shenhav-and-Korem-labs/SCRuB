all: _site.yml
	Rscript -e "rmarkdown::render('index.Rmd')"
	Rscript -e "rmarkdown::render('tutorial.Rmd')"
	Rscript -e "rmarkdown::render('prescrub_intro.Rmd')"
	Rscript -e "rmarkdown::render('Plasma-Tutorial-QIIME2-CLI.Rmd')"
	Rscript -e "rmarkdown::render('qiime2-intro.Rmd')"
	Rscript -e "rmarkdown::render('R-inputs.Rmd')"
