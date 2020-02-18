test_deriv: test_deriv.cpp 
	g++ $^ -o $@

clean:
	rm -f test_deriv
