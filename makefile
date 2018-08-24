CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas

$(OBJ):
	mkdir $(OBJ)

all: \
	permissions \
	BOOM/libBOOM.a \
	swift

BOOM/libBOOM.a:
	cd BOOM ; make all ; cd -

permissions:
	chmod a+x *.pl ; chmod a+x perl/*.pl

#---------------------------------------------------------
$(OBJ)/swift.o:\
		swift.C
	$(CC) $(CFLAGS) -o $(OBJ)/swift.o -c \
		swift.C
#---------------------------------------------------------
swift: \
		$(OBJ) \
		$(OBJ)/Replicate.o \
		$(OBJ)/Replicates.o \
		$(OBJ)/SwiftSample.o \
		$(OBJ)/swift.o
	$(CC) $(LDFLAGS) -o swift \
		$(OBJ)/Replicate.o \
		$(OBJ)/Replicates.o \
		$(OBJ)/SwiftSample.o \
		$(OBJ)/swift.o \
		$(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/Replicate.o:\
		Replicate.C\
		Replicate.H
	$(CC) $(CFLAGS) -o $(OBJ)/Replicate.o -c \
		Replicate.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/Replicates.o:\
		Replicates.C\
		Replicates.H
	$(CC) $(CFLAGS) -o $(OBJ)/Replicates.o -c \
		Replicates.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/SwiftSample.o:\
		SwiftSample.C\
		SwiftSample.H
	$(CC) $(CFLAGS) -o $(OBJ)/SwiftSample.o -c \
		SwiftSample.C
#---------------------------------------------------------
