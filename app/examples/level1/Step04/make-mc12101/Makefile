ROOT   = ../../../../..
TARGET = main.abs

NMCC     = y
NEURO   ?= $(HOME)/nmc4/Bin_old
MC12101 ?= $(HOME)/nmc4/Bin

.PHONY: 1
ifeq ($(OS),Windows_NT)
SHELL=cmd
else
RUNDIR ?= /opt/share/mc12101-proton
PATH := $(PATH):$(NEURO)/bin
endif

LIBS   = mc12101lib_nm.lib libc4f.lib libc4.lib nmblas-nmc4f.lib
CFG    = mc12101brd-nmc0.cfg


$(TARGET): ../main.c ../clock.s Makefile mc12101brd-nmc0.cfg
	asm ../clock.s -nmc4
	nmcc ../main.c -Sc -nmc4 -Tc99 -I$(NEURO)/include -I$(ROOT)/include
	linker clock.elf main.elf -o$(TARGET) -c$(CFG) $(LIBS) -l$(NEURO)/lib -l$(MC12101)/lib -l$(ROOT)/lib -mout.map
	

#	nmcpp -6405 ../main.cpp -OPT0 -debug 
#	asm -nmc4 ../main.asm -I"%NEURO%\include"
#	asm -nmc4 ../add.asm -I"%NEURO%\include"
#	linker main.elf add.elf mc12101lib_nm.lib libc05.lib -l$(MC12101)\lib -l$(NEURO)\lib -cmc12101brd.cfg -m -d0




clean: 
	rm $(TARGET) *.elf *.map


# don't use 'run: $(TARGET)' for remote 	
run: 
	mc12101run $(TARGET) -t10000
	
run-proton: $(TARGET)
	$(RUNDIR)/run

reset:
	pdu.exe 10.7.9.120 user User98 7 3 6

reset-proton:
	$(RUNDIR)/reset12101
