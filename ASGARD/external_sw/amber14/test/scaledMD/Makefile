# Makefile for scaled molecular dynamics tests

all: test

pme:
	@echo "Testing scaledMD with PME"
	@./Run.scaledMD $(OPT)

test: pme 

pmemd: 
	-$(MAKE) test OPT=pmemd
