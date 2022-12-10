all: mainOpenmp.cpp hrtime.h
	rm -f test/test_sequential.txt
	rm -f test/test_OpenMP.txt
	g++ -o sequential -fopenmp sequential.cpp
	g++ -o main -fopenmp mainOpenmp.cpp
clean:
	rm main
	rm sequential