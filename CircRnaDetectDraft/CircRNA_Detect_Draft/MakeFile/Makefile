SharedLibraryPath=../../ShareLibrary/
#BAMTOOLS=$(SharedLibraryPath)bamtools/include
#BAMTOOLS_LD=$(SharedLibraryPath)bamtools/lib

#RootPath = ../ScaffoldGapFillerTest/
#NeedleManPath = ../Depends/needleman_wunsch/
#NormDependsPath = ../Depends/
#INC=-I$(SharedLibraryPath)
#INC_PARAMS=$(foreach d, $(INC), -I$d)
#Notice: "-c" will cause the error

CC = g++
INC_DIR = $(SharedLibraryPath) #Identify include path
#CFLAGS = -O3 -Wall -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
CFLAGS = -O3 -Wall -static
LIBS += -lz
LIBS += -lpthread
#INCLUDEPATH += tr1/memory

OBJ = ./obj/

all: CircRnaDetectDraft

clsbasealgorithm.o : $(SharedLibraryPath)clsbasealgorithm.cpp $(SharedLibraryPath)clsbasealgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbasealgorithm.cpp

clsfastareader.o: $(SharedLibraryPath)clsfastareader.cpp $(SharedLibraryPath)clsfastareader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastareader.cpp

clsfastqreader.o: $(SharedLibraryPath)clsfastqreader.cpp $(SharedLibraryPath)clsfastqreader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastqreader.cpp

#clsbwa.o : $(SharedLibraryPath)clsbwa.cpp $(SharedLibraryPath)clsbwa.h
#	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbwa.cpp

clsgtfparse.o: $(SharedLibraryPath)clsgtfparse.cpp $(SharedLibraryPath)clsgtfparse.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsgtfparse.cpp

clsreadconfigini.o: $(SharedLibraryPath)clsreadconfigini.cpp $(SharedLibraryPath)clsreadconfigini.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsreadconfigini.cpp

clskmeralgorithm.o: $(SharedLibraryPath)clskmeralgorithm.cpp $(SharedLibraryPath)clskmeralgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clskmeralgorithm.cpp

clsconfig.o: ../CircRnaDetectDraft/clsconfig.cpp ../CircRnaDetectDraft/clsconfig.h
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/clsconfig.cpp

clsfindcandidate.o: ../CircRnaDetectDraft/clsfindcandidate.cpp ../CircRnaDetectDraft/clsfindcandidate.h
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/clsfindcandidate.cpp

clskmertable.o: ../CircRnaDetectDraft/clskmertable.cpp ../CircRnaDetectDraft/clskmertable.h
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/clskmertable.cpp

clsresultcomparison.o: ../CircRnaDetectDraft/clsresultcomparison.cpp ../CircRnaDetectDraft/clsresultcomparison.h
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/clsresultcomparison.cpp

clstroubleshoot.o: ../CircRnaDetectDraft/clstroubleshoot.cpp ../CircRnaDetectDraft/clstroubleshoot.h
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/clstroubleshoot.cpp

main.o: ../CircRnaDetectDraft/main.cpp 
	$(CC) $(CFLAGS) -c ../CircRnaDetectDraft/main.cpp

#CircRnaDetectDraft: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsfindcandidate.o clskmertable.o clsresultcomparison.o clstroubleshoot.o main.o 
#	$(CC) $(CFLAGS) -o CircRnaDetectDraft clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsfindcandidate.o clskmertable.o clsresultcomparison.o clstroubleshoot.o main.o \
#	-lbamtools -lz

CircRnaDetectDraft: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsfindcandidate.o clskmertable.o clsresultcomparison.o clstroubleshoot.o main.o 
	$(CC) $(CFLAGS) -o CircRnaDetectDraft clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsfindcandidate.o clskmertable.o clsresultcomparison.o clstroubleshoot.o main.o -lz -lpthread


RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


