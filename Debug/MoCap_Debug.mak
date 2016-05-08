SHELL = cmd.exe

#
# ZDS II Make File - MoCap project, Debug configuration
#
# Generated by: ZDS II - Z8 Encore! Family 4.11.0 (Build 08052301)
#   IDE component: c:4.11:08050501
#   Install Path: C:\Program Files\ZiLOG\ZDSII_Z8Encore!_4.11.0\
#

RM = del

ZDS = C:\PROGRA~1\ZiLOG\ZDSII_~1.0
BIN = $(ZDS)\bin
# ZDS include base directory
INCLUDE = C:\PROGRA~1\ZiLOG\ZDSII_~1.0\include
# intermediate files directory
WORKDIR = Z:\SOFTWA~1\Embedded\MoCap\Debug

CC = @$(BIN)\eZ8cc
AS = @$(BIN)\eZ8asm
LD = @$(BIN)\eZ8link
AR = @$(BIN)\eZ8lib

CFLAGS =  \
-bfpack:tight -fastcall -const:RAM -define:_Z8F6423  \
-define:_Z8ENCORE_XP_64XX_SERIES -define:_Z8ENCORE_F642X  \
-genprintf -NOkeepasm -keeplst -NOlist -NOlistinc -model:S  \
-NOoptlink -promote -regvar:8 -reduceopt  \
-stdinc:"Z:\SOFTWA~1\Embedded\MoCap;$(INCLUDE)\std;$(INCLUDE)\zilog;$(INCLUDE)\zilog\Z8Encore_F642X"  \
-usrinc:"Z:\SOFTWA~1\Embedded\MoCap;" -debug -NOrevaa  \
-cpu:Z8F6423  \
-asmsw:" -cpu:Z8F6423 -define:_Z8F6423=1 -define:_Z8ENCORE_XP_64XX_SERIES=1 -define:_Z8ENCORE_F642X=1 -include:Z:\SOFTWA~1\Embedded\MoCap;$(INCLUDE)\std;$(INCLUDE)\zilog;$(INCLUDE)\zilog\Z8Encore_F642X -NOrevaa"

ASFLAGS =  \
-define:_Z8F6423=1 -define:_Z8ENCORE_XP_64XX_SERIES=1  \
-define:_Z8ENCORE_F642X=1  \
-include:"Z:\SOFTWA~1\Embedded\MoCap;$(INCLUDE)\std;$(INCLUDE)\zilog;$(INCLUDE)\zilog\Z8Encore_F642X"  \
-list -NOlistmac -name -pagelen:56 -pagewidth:80 -quiet -sdiopt  \
-warn -debug -NOigcase -NOrevaa -cpu:Z8F6423

LDFLAGS = @.\MoCap_Debug.linkcmd
OUTDIR = Z:\SOFTWA~1\Embedded\MoCap\Debug

build: MoCap

buildall: clean MoCap

relink: deltarget MoCap

deltarget: 
	@if exist $(WORKDIR)\MoCap.lod  \
            $(RM) $(WORKDIR)\MoCap.lod
	@if exist $(WORKDIR)\MoCap.hex  \
            $(RM) $(WORKDIR)\MoCap.hex
	@if exist $(WORKDIR)\MoCap.map  \
            $(RM) $(WORKDIR)\MoCap.map

clean: 
	@if exist $(WORKDIR)\MoCap.lod  \
            $(RM) $(WORKDIR)\MoCap.lod
	@if exist $(WORKDIR)\MoCap.hex  \
            $(RM) $(WORKDIR)\MoCap.hex
	@if exist $(WORKDIR)\MoCap.map  \
            $(RM) $(WORKDIR)\MoCap.map
	@if exist $(WORKDIR)\main.obj  \
            $(RM) $(WORKDIR)\main.obj
	@if exist $(WORKDIR)\spi.obj  \
            $(RM) $(WORKDIR)\spi.obj

# pre-4.11.0 compatibility
rebuildall: buildall 

LIBS = 

OBJS =  \
            $(WORKDIR)\main.obj  \
            $(WORKDIR)\spi.obj

MoCap: $(OBJS)
	 $(LD) $(LDFLAGS)

$(WORKDIR)\main.obj :  \
            Z:\SOFTWA~1\Embedded\MoCap\src\main.c  \
            $(INCLUDE)\std\STDARG.H  \
            $(INCLUDE)\zilog\FORMAT.H  \
            $(INCLUDE)\zilog\defines.h  \
            $(INCLUDE)\zilog\dmadefs.h  \
            $(INCLUDE)\zilog\ez8.h  \
            $(INCLUDE)\zilog\gpio.h  \
            $(INCLUDE)\zilog\uart.h  \
            $(INCLUDE)\zilog\uartdefs.h  \
            $(INCLUDE)\std\STDARG.H  \
            $(INCLUDE)\std\STDIO.H  \
            $(INCLUDE)\zilog\FORMAT.H  \
            $(INCLUDE)\zilog\SIO.H  \
            $(INCLUDE)\zilog\Zconst.h  \
            $(INCLUDE)\zilog\defines.h  \
            $(INCLUDE)\zilog\ez8.h  \
            Z:\SOFTWA~1\Embedded\MoCap\src\spi.h
	 $(CC) $(CFLAGS) Z:\SOFTWA~1\Embedded\MoCap\src\main.c

$(WORKDIR)\spi.obj :  \
            Z:\SOFTWA~1\Embedded\MoCap\src\spi.c  \
            $(INCLUDE)\zilog\defines.h  \
            $(INCLUDE)\zilog\dmadefs.h  \
            $(INCLUDE)\zilog\gpio.h  \
            $(INCLUDE)\zilog\uart.h  \
            $(INCLUDE)\zilog\uartdefs.h  \
            $(INCLUDE)\zilog\ez8.h  \
            Z:\SOFTWA~1\Embedded\MoCap\src\spi.h
	 $(CC) $(CFLAGS) Z:\SOFTWA~1\Embedded\MoCap\src\spi.c

