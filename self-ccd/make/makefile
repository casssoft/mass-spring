CC=g++
CFLAGS=-include src/vcincludes.h -Iinc -Isrc

all : libselfccd.a

src/%.o : src/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

wildcard := $(wildcard src/*.cpp)
objects := $(patsubst src/%.cpp,src/%.o,$(wildcard src/*.cpp))

libselfccd.a : $(objects)
	ar rvs libselfccd.a $(objects)

$(info $$var is [${objects}])
$(info $$var is [${wildcard}])
