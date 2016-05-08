/////////////////////////////////////////////////////////
/// File: accel.h
///
/// <author>Jeff McCartney</author>
///
/// <version>0.1</version>
///
/// <summary>Defines functions to interact with the
/// accelerometers.  Also, there are definitions for
/// variables that will be needed to support the use of the
/// accelerometers.</summary>
/////////////////////////////////////////////////////////

int g_DataLine = -1;
int g_ClockLine = -1;

/////////////////////////////////////////////////////////
/// Function: setConfig
///
/// <summary>Sets the two lines that are consistant between
/// all accelerometers.  These two lines are the data line
/// and the clock line.</summary>
///
/// <param>dataLine: an integer representing the pin that
/// has been used for the data line between the micro 
/// controller and all the accelerometers.</param>
///
/// <param>clcokLine: an interger representing the pin that
/// has been used for the clcok line between the micro
/// controller and all the accelerometers.</param>
///
/// <returns>Returns 0 for a success. Otherwise the
/// value will be in the error lookup table: <see></see>
/// </returns>
/////////////////////////////////////////////////////////
int setConfig(int dataLine, int clockLine);

/////////////////////////////////////////////////////////
/// Function: getAccelData
///
/// <summary>This function will poll the accelerometer
/// denoted by the passed in ID value.  This function then
/// returns the x,y,z values obtained from the accelerometer.
/// </summary>
///
/// <param>accelID: an integer representination of the pin
/// that triggers the listening of the accelerometer.</param>
///
/// <returns>Returns a string representing the x,y,z values.
/// Each value will be separated by commas.</returns>
/////////////////////////////////////////////////////////
char* getAccelData(int accelID);