src=$(wildcard *.c)
dir=$(notdir $(src))
#dir:=$(notdir ./*.c)
obj=$(dir:%.c=%)
$(obj):$(dir)
	gcc -o $@ $<
all:
	@echo $(dir)
	@echo $(obj)
clean:
	rm $(obj)
