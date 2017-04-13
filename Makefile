exec    = a.out
cc      = g++
opt     = -g
cflags  = -std=c++11 $(opt)

main    = Main.cpp
objects = $(patsubst %.cpp,%.o,$(filter-out $(main), $(wildcard *.cpp)))

.PHONY : all clean

all :	$(objects) 
	@rm -f $(exec)
	@$(MAKE) $(exec)

%.o : %.cpp %.h
	$(cc) $(cflags) -c $<

$(exec) : $(main) $(objects)
	$(cc) $(cflags) $(objects) $< -o $@

clean :
	rm -f $(objects) $(exec)
