simoutbreak.o: simoutbreak.cpp epiworld.hpp
	g++ -std=c++17 -O2 -Wall -pedantic simoutbreak.cpp -o simoutbreak.o

networks: networks.R
	R CMD BATCH --vanilla networks.R networks.Rout &

update:
	rsync -avz ../../world-epi/epiworld.hpp .

sims:
	./simoutbreak.o 1000 ergm && \
		./simoutbreak.o 1000 degseq


.PHONY: networks update sims

