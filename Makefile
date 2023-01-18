#besdtool Makefile

CC := gcc
CPPFLAGS :=
CFLAGS := -g -O0
LDFLAGS :=
LIBS := -lbesd

objs = $(patsubst %.c,%.o,$(wildcard src/*.c))
libobjs = $(patsubst %.c,%.o,$(wildcard lib/*.c))

.PHONY: all
all: besdtool

besdtool: $(objs) libbesd
	@echo besdtool
	$(CC) $(CFLAGS) $(LDFLAGS) $(objs) $(LIBS) -o $@

$(objs): %.o: %.c
	echo $(objs)
	$(CC) $(CPPFLAGS) -c $< -o $@


libbesd: $(libobjs)
	@echo libbesd
	ar rcs libbesd.a $(libobjs)
	

$(libobjs): %.o: %.c
	$(CC) $(CPPFLAGS) -c $< -o $@
