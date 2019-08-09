#
#***********************************************************************
#                                                                      *
#            NACHOS II - A Finite Element Computer Program             *
#                        for Incompressible Flow Problems              *
#                                                                      *
#     Copyright (c) 1986,2019   National Technology & Engineering      *
#                               Solutions of Sandia, LLC (NTESS)       *
#                                                                      *
#                            All rights reserved.                      *
#                                                                      *
#     This software is distributed under the BSD 3-Clause License.     *
#                                                                      *
#***********************************************************************
#
#  makefile for creating nachos ii
#  for linux (7/2008)
#  
#  this is the original version of the code
#
FFLAGS =
CFLAGS = 
CC= icc 
$FC= ifort -cm -r8 -i8
FC= gfortran  -fdefault-real-8 -fdefault-integer-8 

#ACCESS       = /usr/local/eng_sci/struct/current-gcc

LIB_PATH     = -L$(ACCESS)/lib

SUPES_LIB       =  -lsupes
SUPPORT_LIBS    = $(SUPES_LIB) -lsuplib





#
# modules for nachos ii
#
PROGRAM = nachos2.x

OBJS = assmbl.o \
assmkf.o \
assmt.o \
assmu1.o \
assmu2.o \
assmvi.o \
banner.o \
bcksub.o \
bctime.o \
bndry.o \
clsfil.o \
cnnct.o \
convbc.o \
coord3.o \
coord4.o \
coord6.o \
coord8.o \
coord9.o \
deform.o \
diskwr.o \
elemnt.o \
elordr.o \
error.o \
extrap.o \
extrct.o \
fillq.o \
fillt.o \
flux6.o \
flux8.o \
flux9.o \
fluxbc.o \
fluxes.o \
fndpiv.o \
formkf.o \
front.o \
gather.o \
genrd1.o \
genrd2.o \
genrd3.o \
genwrt.o \
grphfl.o \
initlz.o \
invrt3.o \
isrcheq.o \
matrl.o \
maxmin.o \
mesh.o \
midsid.o \
nachos.o \
opnfil.o \
output.o \
point.o \
predct.o \
prevue.o \
prfrnt.o \
prtlmt.o \
ptfunc.o \
quadb.o \
quadc.o \
radbc.o \
radmod.o \
rdffld.o \
rdinpt.o \
relax.o \
sdot.o \
seg1.o \
shape6.o \
shape8.o \
shape9.o \
snorm.o \
solve.o \
spoint.o \
srchl.o \
srchq.o \
store.o \
stream.o \
tnorm.o \
trib.o \
trnsfm.o \
tstep.o \
update.o \
usrbc1.o \
usrbc2.o \
usrbc3.o \
usrbc4.o \
usrbc5.o \
usrbc6.o \
usrcon.o \
usrcp.o \
usrdif.o \
usrex.o \
usrext.o \
usrgrv.o \
usrhtc.o \
usrmsh.o \
usrvhs.o \
usrvis.o \
usrvs.o 



nachos: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(PROGRAM) $(LIB_PATH) $(SUPPORT_LIBS) -lm

clean:
	-rm  -f *.o *~


FFILES = nachos.f \
assmbl.f \
assmkf.f \
assmt.f \
assmu1.f \
assmu2.f \
assmvi.f \
banner.f \
bcksub.f \
bctime.f \
bndry.f \
clsfil.f \
cnnct.f \
convbc.f \
coord3.f \
coord4.f \
coord6.f \
coord8.f \
coord9.f \
deform.f \
diskwr.f \
elemnt.f \
elordr.f \
error.f \
extrap.f \
extrct.f \
fillq.f \
fillt.f \
flux6.f \
flux8.f \
flux9.f \
fluxbc.f \
fluxes.f \
fndpiv.f \
formkf.f \
front.f \
gather.f \
genrd1.f \
genrd2.f \
genrd3.f \
genwrt.f \
grphfl.f \
initlz.f \
invrt3.f \
isrcheq.f \
matrl.f \
maxmin.f \
mesh.f \
midsid.f \
opnfil.f \
output.f \
point.f \
predct.f \
prevue.f \
prfrnt.f \
prtlmt.f \
ptfunc.f \
quadb.f \
quadc.f \
radbc.f \
radmod.f \
rdffld.f \
rdinpt.f \
relax.f \
sdot.f \
seg1.f \
shape6.f \
shape8.f \
shape9.f \
snorm.f \
solve.f \
spoint.f \
srchl.f \
srchq.f \
store.f \
stream.f \
tnorm.f \
trib.f \
trnsfm.f \
tstep.f \
update.f \
usrbc1.f \
usrbc2.f \
usrbc3.f \
usrbc4.f \
usrbc5.f \
usrbc6.f \
usrcon.f \
usrcp.f \
usrdif.f \
usrex.f \
usrext.f \
usrgrv.f \
usrhtc.f \
usrmsh.f \
usrvhs.f \
usrvis.f \
usrvs.f


onefile:
	-cat $(FFILES) > nachos2.for

