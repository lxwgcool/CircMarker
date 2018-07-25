TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += ../../ShareLibrary/
LIBS += -lz
CONFIG += thread

SOURCES += main.cpp \
    ../../ShareLibrary/clsreadconfigini.cpp \
    ../../ShareLibrary/clsfastareader.cpp \
    ../../ShareLibrary/clsbwa.cpp \
    ../../ShareLibrary/clsfastqreader.cpp \
    ../../ShareLibrary/clsbasealgorithm.cpp \
    ../../ShareLibrary/clskmeralgorithm.cpp \
    ../../ShareLibrary/clsgtfparse.cpp \
    clskmertable.cpp \
    clsconfig.cpp \
    clsfindcandidate.cpp \
    clsresultcomparison.cpp \
    clstroubleshoot.cpp

OTHER_FILES += \
    ../Readme/Architecture.txt \
    ../Readme/log_05_11_2016.txt \
    ../../ShareLibrary/fasta_to_fastq.pl \
    ../Readme/RealDataTesting.txt \
    ../Readme/Different_Gene_Version_Description.txt

HEADERS += \
    ../../ShareLibrary/clsfastareader.h \
    ../../ShareLibrary/clsbwa.h \
    ../../ShareLibrary/clsfastqreader.h \
    ../../ShareLibrary/clsbasealgorithm.h \
    ../../ShareLibrary/clskmeralgorithm.h \
    ../../ShareLibrary/clsgtfparse.h \
    clskmertable.h \
    clsconfig.h \
    clsfindcandidate.h \
    ../../ShareLibrary/clsreadconfigini.h \
    clsresultcomparison.h \
    clstroubleshoot.h


unix:!macx: LIBS += -L$$PWD/../../ShareLibrary/bamtools/lib/ -lbamtools

INCLUDEPATH += $$PWD/../../ShareLibrary/bamtools/include
DEPENDPATH += $$PWD/../../ShareLibrary/bamtools/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../ShareLibrary/bamtools/lib/libbamtools.a
