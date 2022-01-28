all: run

# implicit rules
# variables
# $^ : todas las depencias
# $< : Solamente la primera dependencia
# $@ : target

CXX_FLAGS = -O3

%.x : %.cpp
	g++ $(CXX_FLAGS) $< -o $@

run:	Temperatura_Debye.x
	./$<

clean:
	rm -f *.x
