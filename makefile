sol1:
	g++ -O3 --std=c++11 solution-step1.c -o build/sol1
sol2:
	g++ -O3 --std=c++11 solution-step2.c -o build/sol2
sol3:
	g++ -O3 --std=c++11 solution-step3.c -o build/sol3
sol4:
	g++ -O3 --std=c++11 solution-step4.c -o build/sol4
sol5: # pdf
	pdflatex sol4.tex -output-directory latex/
sol6:
	g++ -O3 --std=c++11 solution-step5.c -o build/sol6
