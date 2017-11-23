
.PHONY : dump

CFG = emu.cfg
#CFG =emu1.cfg
TARGET = a.out
$(TARGET): Makefile $(CFG)
	nmc-g++ ../accmul.s ../clock.s sim_utility.s -x c ../main.c -x none  ~/nmc4/Bin/lib/gcc/nmc/4.8.3/../../../../nmc/lib/crt0.o ~/nmc4/Bin/lib/gcc/nmc/4.8.3/../../../../nmc/lib/crtn.o -nostartfiles --entry=super_start -Xlinker  -T $(CFG)
	

run: 
	mc12101run $(TARGET)

run-proton:
	cp * /opt/share/mc12101-proton
	
#mc12101run-proton 
		
emurun: $(TARGET)
	qemu-ppdl -g $(PORT) $(TARGET) &
	nm-gdb 	-ex='target remote localhost:$(PORT)' -q \
		-ex='set disassemble-next-line on' \
		-ex='set confirm off' \
		-ex='load $(TARGET)' \
		-ex='file $(TARGET)' \
		-ex='b*(start+6)' \
		-ex='cont' \
		-ex='p $$gr7' \
		-ex='quit'

verilog: $(TARGET)
	../../bin/nmelf2txt a.out H_E_X.hex tt.hex -f32
	cp H_E_X.hex ~/verilog
	rm ./tt.hex

dump: $(TARGET)
	nmc-objdump -x $(TARGET) > dump.txt

clean:
	rm a.out *~ *.hex *.o dump dump.txt $(TARGET)
