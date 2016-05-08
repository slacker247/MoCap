/////////////////////////////////////////////////////////
/// File: accel.c
///
/// <author>Jeff McCartney</author>
///
/// <version>0.1</version>
///
/// <summary>The accel.h file's function definitions.</summary>
/////////////////////////////////////////////////////////

#include "accel.h"

int setConfig(int dataLine, int clockLine)
{
	int status = 0;
	g_DataLine = dataLine;
	g_ClockLine = clockLine;
	return status;
}

char[] getAccelData(int accelID)
{
	char[] values = "xxx.xx,yyy.yy,zzz.zz";
	// TODO : set accelID pin to low
	// TODO : get read the accel data from the g_DataLine
	// TODO : set accelID pin to high
	// TODO : depending on how the data is returned from
	// the accel format it to fit this "xxx.xx,yyy.yy,zzz.zz".
	return values;
}