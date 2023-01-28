all:
	g++ -std=c++17 -Wall -Wno-unused-function -Wno-sign-compare -O2 splender.cpp -o splender
#g++ -std=c++17 -ggdb -fsanitize=address   -Wall -Wno-unused-function -Wno-sign-compare -O2 splender.cpp -o splender
clean:
	rm splender
