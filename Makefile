bin/bim-innner-join:
	g++ src/bim_inner_join.cpp -o bin/bim-inner-join -std=c++11 -O3

clean:
	rm -f bin/bim-inner-join
