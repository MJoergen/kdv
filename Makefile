kdv: kdv.cpp deriv.cpp integrator.cpp kdvintegrator.cpp
	g++ -Wall -Wextra -O3 $^ -o $@

test_deriv: test_deriv.cpp deriv.cpp
	g++ -Wall -Wextra -O3 $^ -o $@

clean:
	rm -f kdv
	rm -f test_deriv
