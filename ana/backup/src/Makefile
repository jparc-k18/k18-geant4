#
#	Makefile for Sks Analyzers
#
#

CC	= gcc
C++	= g++

CLHEPDIR        = /common/CLHEP/2.0.4.2
#CLHEPDIR        = /sw/packages/CLHEP/2.1.1.0
#
CLHEPINCDIR     = $(CLHEPDIR)/include
CLHEPLIBDIR     = $(CLHEPDIR)/lib
#
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIBS	= -L$(CLHEPLIBDIR) -lCLHEP

RTINCLUDE       = `root-config --cflags`
RTLIBS          = `root-config --libs`

### HD DAQ Unpacker ####
unpacker_install      := $(HOME)/work/K1.8/K1.8Analyzer/unpacker
unpacker_include_dir  := $(unpacker_install)/include
unpacker_cflags       := -I$(unpacker_include_dir)
unpacker_lib_dir      := $(unpacker_install)/lib
unpacker_libs         := \
 -Wl,-rpath,$(unpacker_lib_dir) -L$(unpacker_lib_dir) -lHDDAQUnpacker


CPPFLAGS	= $(DEFINEFLAGS) $(DEBUGFLAGS) -I.	\
		$(CLHEPINCLUDE) $(RTINCLUDE) 	\
		$(unpacker_cflags)

CFLAGS		= -g -O2


LIBS		= $(CLHEPLIBS) $(RTLIBS)

RM	= rm -f
TOUCH	= touch

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

.cc.o:
	$(C++) $(CFLAGS) $(CPPFLAGS) -c $<

.cpp.o:
	$(C++) $(CFLAGS) $(CPPFLAGS) -c $<

all:
	make lib; make programs


LIBOBJS	= DCHit.o DCLocalTrack.o DCConfMan.o TrackHit.o		\
	FieldMan.o SksFieldMap.o RungeKuttaUtilities.o		\
	FieldElements.o DCGeomRecord.o SksTrack.o MathTools.o   \
	Kinematics.o FiberHit.o FiberTrack.o DCAnalyzer.o       \
	DCTrackSearch.o CylindricalAna.o EnergyCorrection.o     \
	EvDisp.o Minuit.o K18Track.o K18TrackFCN.o K18TrackD2UFCN.o \
	K18TransMatrix.o

LIBARC	= libSksAna.a

SKSLIB	= -L. -lSksAna

BINDIR 	= ../bin

PROGRAMS = ../bin/analysSksTracking ../bin/analysSksTrackingForYd	\
	   ../bin/analysLocalTrackingIn ../bin/analysLocalTrackingOut 	\
	   ../bin/analysDCConfMan ../bin/analysFieldMan \
	   ../bin/analysCylindricalAna  \
	   ../bin/analysSksTrackingForSigmaPlusP \
	   ../bin/analysCylindricalAnaSigmaPlus \
	   ../bin/analysCylindricalEvDisp \
	   ../bin/StsTracking 

lib: $(LIBARC)

$(LIBARC): $(LIBOBJS)
	ar r $(LIBARC) $(LIBOBJS) 

programs: $(PROGRAMS) $(LIBARC)

clean:
	$(RM) *~ *.o $(LIBOBJS) $(LIBARC) $(PROGRAMS)

# Programs

$(BINDIR)/analysSksTrackingForYd: UserSksTrackingForYd.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserSksTrackingForYd.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysSksTrackingForSigmaPlusP: UserSksTrackingForSigmaPlusP.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserSksTrackingForSigmaPlusP.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysSksTracking: UserSksTracking.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserSksTracking.o $(SKSLIB) $(LIBS)
$(BINDIR)/StsTracking: UserSksTracking.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserSksTracking.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysLocalTrackingIn: UserLocalTrackingIn.o $(LIBARC) 
	$(C++) -o $@ $(CFLAGS) UserLocalTrackingIn.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysLocalTrackingOut: UserLocalTrackingOut.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserLocalTrackingOut.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysDCConfMan: UserDCConfMan.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserDCConfMan.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysFieldMan: UserFieldMan.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserFieldMan.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysCylindricalAna: UserCylindricalAna.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserCylindricalAna.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysCylindricalAnaSigmaPlus: UserCylindricalAnaSigmaPlus.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserCylindricalAnaSigmaPlus.o $(SKSLIB) $(LIBS)

$(BINDIR)/analysCylindricalEvDisp: UserCylindricalEvDisp.o $(LIBARC)
	$(C++) -o $@ $(CFLAGS) UserCylindricalEvDisp.o $(SKSLIB) $(LIBS)

UserDCConfMan.o: UserDCConfMan.cc DCConfMan.hh ThreeVector.hh

UserLocalTrackingIn.o: UserLocalTrackingIn.cc DCLocalTrack.hh DCHit.hh	\
		DCConfMan.hh RootHelper.hh

UserLocalTrackingOut.o: UserLocalTrackingOut.cc DCLocalTrack.hh 	\
		DCHit.hh DCConfMan.hh RootHelper.hh

UserSksTracking.o: UserSksTracking.cc FieldMan.hh DCLocalTrack.hh	\
		DCHit.hh DCConfMan.hh SksTrack.hh TrackHit.hh RootHelper.hh

UserSksTrackingForYd.o: UserSksTrackingForYd.cc FieldMan.hh DCLocalTrack.hh \
		DCHit.hh DCConfMan.hh SksTrack.hh TrackHit.hh RootHelper.hh

UserSksTrackingForSigmaPlusP.o: UserSksTrackingForSigmaPlusP.cc FieldMan.hh \
		DCLocalTrack.hh DCHit.hh DCConfMan.hh SksTrack.hh \
		TrackHit.hh RootHelper.hh

UserCylindricalAna.o: UserCylindricalAna.cc FieldMan.hh DCLocalTrack.hh \
		DCHit.hh DCConfMan.hh SksTrack.hh TrackHit.hh RootHelper.hh

UserCylindricalAnaSigmaPlus.o: UserCylindricalAnaSigmaPlus.cc FieldMan.hh \
		DCLocalTrack.hh \
		DCHit.hh DCConfMan.hh SksTrack.hh TrackHit.hh RootHelper.hh

UserCylindricalEvDisp.o: UserCylindricalEvDisp.cc FieldMan.hh DCLocalTrack.hh \
		DCHit.hh DCConfMan.hh SksTrack.hh TrackHit.hh RootHelper.hh  \
		EvDisp.hh

UserFieldMan.o: UserFieldMan.cc FieldMan.hh ThreeVector.hh

# Object Files

DCHit.o: DCHit.cc DCHit.hh

DCLocalTrack.o: DCLocalTrack.cc DCLocalTrack.hh DCHit.hh DCConfMan.hh	\
		MathTools.hh

DCConfMan.o: DCConfMan.cc DCConfMan.hh DCGeomRecord.hh TemplateLib.hh

TrackHit.o: TrackHit.cc TrackHit.hh DCConfMan.hh

FieldMan.o: FieldMan.cc FieldMan.hh SksFieldMap.hh FieldElements.hh

SksFieldMap.o: SksFieldMap.cc SksFieldMap.hh

RungeKuttaUtilities.o: RungeKuttaUtilities.cc RungeKuttaUtilities.hh	\
		FieldMan.hh DCConfMan.hh DCGeomRecord.hh

FieldElements.o: FieldElements.cc FieldElements.hh

DCGeomRecord.o: DCGeomRecord.cc DCGeomRecord.hh

SksTrack.o: SksTrack.cc SksTrack.hh DCLocalTrack.hh TrackHit.hh		\
		DCConfMan.hh MathTools.hh

MathTools.o: MathTools.cc MathTools.hh

Kinematics.o: Kinematics.cc Kinematics.hh

FiberHit.o: FiberHit.cc FiberHit.hh

FiberTrack.o: FiberTrack.cc FiberTrack.hh

DCAnalyzer.o: DCAnalyzer.cc DCAnalyzer.hh DetectorId.hh DCHit.hh  \
	      DCLocalTrack.hh SksTrack.hh TemplateLib.hh DCTrackSearch.hh \
	      DCConfMan.hh myfunc.hh

DCTrackSearch.o: DCTrackSearch.cc DCTrackSearch.hh DCParameters.hh \
	         DCLocalTrack.hh TemplateLib.hh DetectorId.hh DCConfMan.hh

CylindricalAna.o: CylindricalAna.cc CylindricalAna.hh YdScatDetectorConst.hh \
	          RadDeg.hh

EnergyCorrection.o: EnergyCorrection.cc EnergyCorrection.hh  \
		YdScatDetectorConst.hh 

EvDisp.o: EvDisp.cc EvDisp.hh  YdScatDetectorConst.hh 

Minuit.o: Minuit.cc Minuit.hh 

K18Track.o: K18Track.cc DCLocalTrack.hh K18TransMatrix.hh DCConfMan.hh \
	DetectorId.hh TrackHit.hh TemplateLib.hh K18TrackFCN.hh \
	K18TrackD2UFCN.hh Minuit.hh K18Parameters.hh DCAnalyzer.hh

K18TrackFCN.o: K18TrackFCN.cc K18TrackFCN.hh K18Track.hh K18TransMatrix.hh \
	TrackHit.hh DCConfMan.hh

K18TrackD2UFCN.o: K18TrackD2UFCN.cc K18TrackD2UFCN.hh K18Track.hh \
	K18TransMatrix.hh TrackHit.hh DCConfMan.hh

K18TransMatrix.o: K18TransMatrix.cc K18TransMatrix.hh

# Include Files

DCHit.hh:

DCLocalTrack.hh:

DCConfMan.hh: ThreeVector.hh
	$(TOUCH) $@

DCGeomRecord.hh: ThreeVector.hh
	$(TOUCH) $@

TrackHit.hh: DCHit.hh ThreeVector.hh
	$(TOUCH) $@

FieldMan.hh: ThreeVector.hh
	$(TOUCH) $@

SksFieldMap.hh: 

RungeKuttaUtilities.hh: ThreeVector.hh
	$(TOUCH) $@

FieldElements.hh: DCGeomRecord.hh ThreeVector.hh
	$(TOUCH) $@

SksTrack.hh: RungeKuttaUtilities.hh ThreeVector.hh
	$(TOUCH) $@


MathTools.hh:

ThreeVector.hh:

TemplateLib.hh:

Kinematics.hh: ThreeVector.hh
	$(TOUCH) $@

FiberHit.hh: 

FiberTrack.hh: 

DCAnalyzer.hh: DetectorId.hh ThreeVector.hh

DCParameters.hh: 

DCTrackSearch.hh: DCAnalyzer.hh ThreeVector.hh

DetectorId.hh: 

DadDeg.hh:

CylindricalAna.hh: ThreeVector.hh

EnergyCorrection.hh: 

EvDisp.hh: 

Minuit.hh:

K18Track.hh: ThreeVector.hh

K18TrackFCN.hh: Minuit.hh

K18TrackD2UFCN.hh: Minuit.hh

K18TransMatrix.hh: