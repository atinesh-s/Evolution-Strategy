all: main
main: main.o ESSolver.o
	g++ obj/main.o obj/ESSolver.o -o main -pthread
ESSolver.o:
	g++ -c src/ESSolver.cpp -o obj/ESSolver.o -std=c++11 
main.o:
	g++ -c src/main.cpp -o obj/main.o -std=c++11 
