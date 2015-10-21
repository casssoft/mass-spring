CC=g++
CFLAGS=-include src/vcincludes.h -Iinc -Isrc

all : libselfccd.a sample/a.out

src/%.o : src/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

libselfccd.a : src/*.o
	ar rvs libselfccd.a src/*.o

sample/a.out : libselfccd.a
	g++ -include ../src/vcincludes.h *.cpp ../libselfccd.a -I../inc -o sample/a.out

