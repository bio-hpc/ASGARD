# Makefile for accelerated molecular dynamics tests

all: test

pme:
	@echo "Testing AMD with PME"
	@-cd dhfr_pme && ./Run.pme.amd1 $(OPT)
	@-cd dhfr_pme && ./Run.pme.amd2 $(OPT)
	@-cd dhfr_pme && ./Run.pme.amd3 $(OPT)

ips:
	@echo "Testing AMD with IPS"
	@-cd gact_ips && ./Run.ips.amd1 $(OPT)
	@-cd gact_ips && ./Run.ips.amd2 $(OPT)
	@-cd gact_ips && ./Run.ips.amd3 $(OPT)

gb:
	@echo "Testing AMD with GB"
	@-cd rna_gb && ./Run.gb.amd1 $(OPT)
	@-cd rna_gb && ./Run.gb.amd2 $(OPT)
	@-cd rna_gb && ./Run.gb.amd3 $(OPT)

test: pme ips gb

pmemd: 
	-$(MAKE) test OPT=pmemd
