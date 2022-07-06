simoutbreak.o: simoutbreak.cpp epiworld.hpp
	g++ -std=c++17 -g -O2 -Wall -pedantic simoutbreak.cpp -o simoutbreak.o

networks: networks.R
	R CMD BATCH --vanilla networks.R networks.Rout &

update:
	rsync -avz ../../world-epi/epiworld.hpp .

Size=1454

sims: simoutbreak.o
	./simoutbreak.o ergm 1000 ${Size} && \
		./simoutbreak.o degseq 1000 ${Size} && \
		./simoutbreak.o ergm+degseq 1000 ${Size} && \
		./simoutbreak.o original 1000 ${Size} && \
		./simoutbreak.o permute 1000 ${Size} && \
		$(MAKE) -B abm.pdf &

abm.pdf: abm.Rmd simoutbreak-params.txt
	Rscript --vanilla -e 'rmarkdown::render("abm.Rmd")' && evince abm.pdf &

.PHONY: networks update sims

