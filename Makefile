all:
	g++ -std=c++17 -D SPLENDER_DEBUG -Wall -Wno-unused-function -Wno-sign-compare -O2 splender.cpp -o splender
