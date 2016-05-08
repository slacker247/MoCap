// SerialPort.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <windows.h>
#include <conio.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char* argv[])
{
    ofstream myFile;
    char c = '0';
    char INBUFFER[500];
    char OUTBUFFER[20];
    DWORD        bytes_read    = 0;

    // Number of bytes read from port
    DWORD        bytes_written = 0;
    // Number of bytes written to the port
    HANDLE       comport      = NULL;
    // Handle COM port
    int   bStatus;
    DCB          comSettings;
    // Contains various port settings
    COMMTIMEOUTS CommTimeouts;

    LPCWSTR port = _T("\\\\.\\COM4");

    strcpy(&OUTBUFFER[0], "The quick brown fox jumped over the lazy dog. \n\r\0");
    // Open COM port
    if ((comport =
        CreateFile(port,     // open com5:
        GENERIC_READ | GENERIC_WRITE, // for reading and writing
        0,                            // exclusive access
        NULL,                         // no security attributes
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL)) == INVALID_HANDLE_VALUE)
    {
        // error processing code goes here
    }
    // Set timeouts in milliseconds
    CommTimeouts.ReadIntervalTimeout         = 0;
    CommTimeouts.ReadTotalTimeoutMultiplier  = 0;
    CommTimeouts.ReadTotalTimeoutConstant    = 100;
    CommTimeouts.WriteTotalTimeoutMultiplier = 0;
    CommTimeouts.WriteTotalTimeoutConstant   = 100;
    bStatus = SetCommTimeouts(comport,&CommTimeouts);
    if (bStatus != 0)
    {
        // error processing code goes here
    }
    // Set Port parameters.
    // Make a call to GetCommState() first in order to fill
    // the comSettings structure with all the necessary values.
    // Then change the ones you want and call SetCommState().
    GetCommState(comport, &comSettings);
    comSettings.BaudRate = 115200;
    comSettings.StopBits = ONESTOPBIT;
    comSettings.ByteSize = 8;
    comSettings.Parity   = NOPARITY;
    comSettings.fParity  = FALSE;
    bStatus = SetCommState(comport, &comSettings);
    if (bStatus == 0)
    {
        // error processing code goes here
    }
    myFile.open ("Test.txt");
    while(!kbhit())
    {
        //bStatus = WriteFile(comport,              // Handle
        //    &OUTBUFFER,      // Outgoing data
        //    48,              // Number of bytes to write
        //    &bytes_written,  // Number of bytes written
        //    NULL);
        //if (bStatus != 0)
        //{
        //    // error processing code here
        //}
        bStatus = ReadFile(comport,   // Handle
            &INBUFFER,            // Incoming data
            1,                  // Number of bytes to read
            &bytes_read,          // Number of bytes read
            NULL);
        if (bStatus != 0)
        {
            // error processing code goes here
        }
        // code to do something with the data goes here
        for(int i = 0; i < bytes_read; i++)
        {
            if(INBUFFER[i] == 'n')
                INBUFFER[i] = '\n';
        }
        INBUFFER[bytes_read] = '\0';
        //printf("%s", INBUFFER);
        myFile << INBUFFER;
        //myFile.flush();
        memset(INBUFFER, 0, 500);
    }
    myFile.close();
    CloseHandle(comport);
    return 0;
}
