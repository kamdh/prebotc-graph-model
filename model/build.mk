# build the compiled code, not much here because it's mostly python

helpers.o: helpers.c helpers.h
	gcc -c -fPIC -o helpers.o helpers.c

libhelpers.a: helpers.o
	ar rcs libhelpers.a helpers.o