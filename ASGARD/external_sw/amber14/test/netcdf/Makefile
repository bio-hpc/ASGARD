# Makefile for Netcdf Restart tests

all: test

mdrstwrite: AlaDipeptide.TIP3P.parm7
	@-cd MD_ncrestart_write && ./runmd.sh $(OPT)

minrstwrite: AlaDipeptide.TIP3P.parm7
	@-cd Min_ncrestart_write && ./runmin.sh $(OPT)

ncrefmd: AlaDipeptide.TIP3P.parm7
	@-cd refMD_ncrestart && ./runmd.sh $(OPT)

ntx5: AlaDipeptide.TIP3P.parm7
	@-cd restartMD.ntx5 && ./runmd.sh $(OPT)

ntx1: AlaDipeptide.TIP3P.parm7
	@-cd restartMD.ntx1 && ./runmd.sh $(OPT)

test: mdrstwrite minrstwrite ncrefmd ntx5 ntx1 

clean:
	-$(MAKE) test OPT=clean

