

#
# C compiler specification
#
CC = cc
#
# C optimization, only one CFLAGS line should be uncommented
# The -Xcpluscomm flags is required for the SGI cc compiler to recognize the C++
# style comments ( // ) used in the code.
#
#CFLAGS = -Xcpluscomm 
#CFLAGS = -Xcpluscomm -O
#CFLAGS = -Xcpluscomm -O2
#CFLAGS = -Xcpluscomm -O3
#CFLAGS = -Xcpluscomm -n32
#CFLAGS = -Xcpluscomm -O -n32 
#CFLAGS = -Xcpluscomm -O2 -n32
#CFLAGS = -Xcpluscomm -O3 -n32
#CFLAGS = -Xcpluscomm -g
CFLAGS = -g

all:
	@echo "To make an executable, type one of the following:"
	@echo " "
	@echo " make toyPlane"
	@echo " "
	@echo "To test one of the executables, type one of the following:"
	@echo " "
	@echo " make testc"
	@echo " "
	@echo "To remove all the generated data, type the following:"
	@echo " "
	@echo " make cleandata"
	@echo " "
	@echo "To remove all the generated data and compiled code, type the following:"
	@echo " "
	@echo " make clean"
	@echo " "


#
# Creation of ToyFDTD from C source code
#
toyPlane: toyPlane.o hsvrgb.o cmap2.o output.o toyPlane.h Makefile
	${CC} ${CFLAGS} -o toyPlane toyPlane.o hsvrgb.o cmap2.o output.o -lm


toyPlane.o: toyPlane.c toyPlane.h Makefile
	${CC} ${CFLAGS} -c -o toyPlane.o toyPlane.c

output.o: output.c toyPlane.h Makefile
	${CC} ${CFLAGS} -c output.c

hsvrgb.o: hsvrgb.c Makefile
	${CC} ${CFLAGS} -c hsvrgb.c

cmap2.o: cmap2.c Makefile
	${CC} ${CFLAGS} -c cmap2.c




testc: toyPlane cleandata
	rm -f c_runLog
	time toyPlane > c_runLog


#
# Type: make clean
# to remove executables, core files, object files, et cetera
#
clean: cleandata
	rm -f toyPlane *.o *~ core
#
# Type: make cleandata
# to remove only data files created by ToyFDTD
#
cleandata:
	rm -f c_*0.mtv
	rm -f c_*1.mtv
	rm -f c_*2.mtv
	rm -f c_*3.mtv
	rm -f c_*4.mtv
	rm -f c_*5.mtv
	rm -f c_*6.mtv
	rm -f c_*7.mtv
	rm -f c_*8.mtv
	rm -f c_*9.mtv
	rm -f c_*0.pnm
	rm -f c_*1.pnm
	rm -f c_*2.pnm
	rm -f c_*3.pnm
	rm -f c_*4.pnm
	rm -f c_*5.pnm
	rm -f c_*6.pnm
	rm -f c_*7.pnm
	rm -f c_*8.pnm
	rm -f c_*9.pnm
	rm -f c_runLog
	rm -f ToyPlaneFDTDc.viz
	rm -f c_*.dat




.PHONY : clean cleandata all depend testc





