all:
	g++ -std=c++17 -D SPLENDER_INTERACTIVE -Wall -Wno-unused-function -O2 splender.cpp -o splender
