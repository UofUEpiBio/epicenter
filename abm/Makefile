simoutbreak.o: simoutbreak.cpp epiworld.hpp
	g++ -std=c++17 -O2 -Wall -pedantic simoutbreak.cpp -o simoutbreak.o

networks: networks.R
	R CMD BATCH --vanilla networks.R networks.Rout &

update:
	rsync -avz ../../world-epi/epiworld.hpp .

sims:
	./simoutbreak.o ergm 1000 981 && \
		./simoutbreak.o degseq 1000 981 && \
		Rscript --vanilla -e 'rmarkdown::render("abm.Rmd")' && \
		evince abm.pdf &


.PHONY: networks update sims

