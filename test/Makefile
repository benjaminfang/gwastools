.PHONY: besd_test
besd_test: besdfile_test.c ../lib/besdfile.c
	gcc besdfile_test.c -c
	gcc ../lib/besdfile.c -c
	gcc besdfile_test.o besdfile.o -o besd_test
