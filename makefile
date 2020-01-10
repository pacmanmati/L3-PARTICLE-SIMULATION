base:
	mkdir -p bin/base
	g++ -O3 --std=c++11 assignment-2019.c -o bin/base/base
sol1:
	mkdir -p bin/sol1
	g++ -O3 --std=c++11 solution-step1.c -o bin/sol1/sol1
sol2:
	mkdir -p bin/sol2
	g++ -O3 --std=c++11 solution-step2.c -o bin/sol2/sol2
sol3:
	mkdir -p bin/sol3
	g++ -O3 --std=c++11 solution-step3.c -o bin/sol3/sol3
sol4:
	mkdir -p bin/sol4
	g++ -O3 --std=c++11 solution-step4.c -o bin/sol4/sol4 -fopenmp
sol5: # pdf
	pdflatex sol4.tex -output-directory latex/
sol6:
	mkdir -p bin/sol6
	g++ -O3 --std=c++11 solution-step5.c -o bin/sol6/sol6
