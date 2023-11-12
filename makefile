
FLAG = -O3 -m64 -Wall

Matching.o: ./Matching.cpp ./Matching.h ./BinaryHeap.h ./Globals.h ./Graph.h 
	g++ $(FLAG) -c ./Matching.cpp -o Matching.o

BinaryHeap.o: ./BinaryHeap.h ./BinaryHeap.cpp ./Globals.h
	g++ $(FLAG) -c ./BinaryHeap.cpp -o BinaryHeap.o

Graph.o: ./Graph.h ./Graph.cpp
	g++ $(FLAG) -c ./Graph.cpp -o Graph.o

TSPLIB_parser.o: TSPLIB_parser.cpp TSPLIB_parser.h ./Graph.h 
	g++ $(FLAG) -c TSPLIB_parser.cpp -o TSPLIB_parser.o

Example.o: Example.cpp MST.h ./Matching.h ./Graph.h Christofides.h TSPLIB_parser.h
	g++ $(FLAG) -c Example.cpp -o Example.o

christofides: Matching.o BinaryHeap.o Graph.o Example.o TSPLIB_parser.o
	g++ $(FLAG) Matching.o BinaryHeap.o Graph.o Example.o TSPLIB_parser.o -o christofides
