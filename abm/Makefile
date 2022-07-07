simoutbreak.o: simoutbreak.cpp epiworld.hpp
	g++ -std=c++17 -g -O2 -Wall -pedantic simoutbreak.cpp -o simoutbreak.o

networks: networks.R
	R CMD BATCH --vanilla networks.R networks.Rout &

update:
	rsync -avz ../../world-epi/epiworld.hpp .

Size=507
NNets=1000

sims: simoutbreak.o
	./simoutbreak.o ergm ${NNets} ${Size} && \
		# ./simoutbreak.o degseq ${NNets} ${Size} && \
		./simoutbreak.o ergm+degseq ${NNets} ${Size} && \
		./simoutbreak.o ergm+attrs ${NNets} ${Size} && \
		./simoutbreak.o original ${NNets} ${Size} && \
		# ./simoutbreak.o permute ${NNets} ${Size} && \
		$(MAKE) -B abm.pdf &

abm.pdf: abm.Rmd simoutbreak-params.txt
	Rscript --vanilla -e 'rmarkdown::render("abm.Rmd")' && evince abm.pdf &

.PHONY: networks update sims

