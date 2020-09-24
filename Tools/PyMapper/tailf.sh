#!/bin/bash

for ((;;)) do
	clear
	tail -n15 log.txt | nl -ba
	sleep 1
done
