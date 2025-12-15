all: main.cpp main.o
	g++ -I/home/utep/gsl/include -c main.cpp
	g++ -L/home/utep/gsl/lib main.o -lgsl -lgslcblas -lm

run: a.out
	./a.out


