kdv: kdv.cpp deriv.cpp
	g++ $^ -o $@

clean:
	rm -f kdv
