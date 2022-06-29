simoutbreak.o: simoutbreak.cpp epiworld.hpp
	g++ -std=c++17 -O2 -Wall -pedantic simoutbreak.cpp -o simoutbreak.o
ergm: networks.R
	R CMD BATCH --vanilla networks.R networks.Rout &
.PHONY: ergm
