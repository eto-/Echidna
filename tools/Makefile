# Postgres libs
ifdef PGLIB
LDFLAGS += -L$(PGLIB) -lpq
else
LDFLAGS += -lpq
endif
ifdef PGINCLUDE
CXXFLAGS += -I$(PGINCLUDE)
else
CXXFLAGS += -I/usr/include/postgresql -I/usr/include/pgsql
endif

# Z libs
LDFLAGS += -lz

all: fix_db

fix_db: fix_db.cc

fake_omon: fake_omon.cc
	g++  fake_omon.cc -o fake_omon -lz

omon_read: omon_read.c
	gcc -Wall -g -o $@ $^
