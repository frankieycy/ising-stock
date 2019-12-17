#!/bin/bash

# clear old files, execute main, generate plot

rm -rf out/* plt/*
clang -o main main.c
if [ $? -eq 0 ]; then
	./main
	python plot_stock.py
fi
