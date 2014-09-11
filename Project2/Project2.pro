TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp
LIBS += -larmadillo -llapack -lblas
QMAKE_CXXFLAGS += -Wall
