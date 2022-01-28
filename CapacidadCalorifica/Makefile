all: run

# implicit rules
# variables
# $^ : todas las depencias
# $< : Solamente la primera dependencia
# $@ : target

CXX_FLAGS = -O3

%.x : %.cpp
	g++ $(CXX_FLAGS) $< -o $@

run:	a.out
	./$<

a.out:	version2_temp_debye.cpp
	g++ -Ofast $<

clean:
	rm -f datos.txt updown.x

clean-all:
	rm -f datos.txt updown.x datos.pdf
