# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -Wno-missing-field-initializers -g
LM = -lm

# optimized compiling
#CX = gcc -O3 -Wall -Wextra -Wno-missing-field-initializers
#LMX = -lm

# using intel C compiler
CX = icc -O2 -Wall -Wremarks
LMX =

LGL = -lGL -lglut

csrc = $(wildcard *.c)
bins = $(patsubst %.c, %, $(csrc))
bins_d = $(patsubst %, %_d, $(bins))
bins_nofftw = $(patsubst %, %_nofftw, $(bins))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LGL)

$(bins_d): %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LGL)

clean:
	rm -rf $(bins) $(bins_d) *.pos *~ a.out *.dat* gmon.out .*.un~ */*~ */.*.un~ MTSEED
	rm -rf r[0-9]*hs
	rstrip.py -Rlv

