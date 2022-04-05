#!/bin/bash
case $1 in
"clean")
	rm -r ~/.R
	;;
*)
	mkdir ~/.R
	touch ~/.R/Makevars
	echo "CXXFLAGS= -O3 -std=c++14 -Wall -march=native" >~/.R/Makevars
	echo "CXX14FLAGS = -O3 -std=c++14 -Wall -march=native" >>~/.R/Makevars
	echo "MAKEFLAGS = -j4" >>~/.R/Makevars
	;;
esac
