# build the compiled code, not much here because it's mostly python

helpers.o: helpers.c helpers.h
	gcc -c -o helpers.o helpers.c