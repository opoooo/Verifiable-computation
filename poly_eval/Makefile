#makefile
test: 
	gcc md2.c lhe.c lhep.c bn_ext.c hash.c poly_eval.c poly_eval_test.c error_hdl.c -I /usr/local/include/flint  -lgmp -lflint -lrelic -lrelic_s -o test
.PHONY : clean
 clean:
	rm -f test
	
