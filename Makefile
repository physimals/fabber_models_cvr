include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_cvr

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_CPROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
  MATLIB = -lnewmat
else 
  UNAME := $(shell uname -s)
  ifeq ($(UNAME), Linux)
    MATLIB = -lopenblas
  endif
  NIFTILIB = -lNewNifti
endif

LIBS = -lnewimage -lmiscmaths -lutils -lprob ${MATLIB} ${NIFTILIB} -lznz -lz -ldl

XFILES = fabber_cvr

# Forward models
OBJS =  fwdmodel_cvr.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all:	${XFILES} libfabbermodels_cvr.a

# models in a library
libfabbermodels_cvr.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_cvr : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
