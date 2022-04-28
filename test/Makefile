test: 
	gcc -g -c test.c
	gcc -g -c ../lib/cut_verylong_ref_alt.c
	gcc -g test.o cut_verylong_ref_alt.o -o test

clean:
	rm -rf *.o
	rm -rf test
