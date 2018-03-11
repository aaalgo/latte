CC = g++
CXXFLAGS += -O3 -g -std=c++11 -I/usr/include/python2.7 -Iaffy/sdk/ -D_INCLUDE_UNISTD_HEADER_=1 -fopenmp
LDFLAGS += -L/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu -L/usr/lib/python2.7/config-x86_64-linux-gnu/ -fopenmp
LDLIBS += -lboost_program_options -lboost_numpy -lboost_python -lpython2.7

AFFY_SRCS = affy/sdk/file/FileIO.cpp \
	    affy/sdk/file/CELFileData.cpp \
	    affy/sdk/file/CDFFileData.cpp 

AFFY_OBJS = $(AFFY_SRCS:.cpp=.o)


all:	pack compare extract_ranks


cdfcel:	cdfcel.o $(AFFY_OBJS)
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


