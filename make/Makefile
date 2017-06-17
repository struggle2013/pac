ALL_FILE=$(notdir ./*.c)
NAME_FILE=$(patsubst %.c,%,$(ALL_FILE))
all:
	@echo $(ALL_FILE)
	@echo $(NAME_FILE)
%:%.c
	gcc -o $@ $<

