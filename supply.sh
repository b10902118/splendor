#!/bin/bash

for n in {1..6}; do
for i in {1..6}; do
	./splender  <<< "$i" | awk 'NF<8' >>"supply_statistic$n" 
done
done
