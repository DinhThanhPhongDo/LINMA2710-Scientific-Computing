SRCS=SparseMatrix.cpp Vector.cpp SparseVector.cpp AbstractVector.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
CPPFLAGS=-Wall -g -std=c++11
LDFLAGS=-Wall -g -std=c++11

runtests: test
	./test

memtest: test
	valgrind ./test

all: test

test: test.o $(OBJS)
	$(CXX) $(LDFLAGS) -o test test.o $(OBJS) $(LDLIBS)

test.o: test.cpp $(SRCS)

SparseMatrix.o: SparseMatrix.cpp SparseMatrix.hpp Vector.hpp SparseVector.hpp AbstractVector.hpp

Vector.o: Vector.cpp Vector.hpp AbstractVector.hpp

SparseVector.o: SparseVector.cpp SparseVector.hpp AbstractVector.hpp

AbstractVector.o: AbstractVector.cpp AbstractVector.hpp

clean:
	$(RM) $(OBJS) *.o *.$(DLEXT)

distclean: clean
	$(RM) test
