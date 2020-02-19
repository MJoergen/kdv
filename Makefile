kdv: kdv.cpp deriv.cpp
	g++ $^ -o $@

test_deriv: test_deriv.cpp deriv.cpp
	g++ $^ -o $@

clean:
	rm -f kdv
	rm -f test_deriv
