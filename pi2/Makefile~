#makefile
test: 
	gcc md2.c lhe.c lhep.c bn_ext.c hash.c pi2.c pi2_test.c error_hdl.c -I /usr/local/include/flint  -lgmp -lflint -lrelic -lrelic_s -o test
.PHONY : clean
 clean:
	rm -f test
	
