/////////////////////////////////////////////////////////
/// File: main.c
///
/// <author>Jeff McCartney</author>
///
/// <version>0.1</version>
///
/// <summary>The main entry point for the application.
/// Pins for the accelerometers are as follows:
/// PF7 - DIO
/// PF6 - CLK
/// PE7 - CS  // High means not listening, low means listening
///
/// PE 0 - 7
/// PF 1 - 7
/// PG 1 - 7
/// </summary>
/////////////////////////////////////////////////////////

#include <ez8.h>
#include <defines.h>
#include <stdio.h>
#include <sio.h>

#include "spi.h"
//#include "accel.h"

//#define SPIPort 0x3C    // PC5=MISO, PC4=MOSI, PC3=SCK, PC2=SS

const int g_NumAccel = 1;

/////////////////////////////////////////////////////////
/// Function: writeData
///
/// <summary>Write the data to the host computer.</summary>
///
/// <returns>Nothing</returns>
/////////////////////////////////////////////////////////
void writeData(char x, char y, char z)
{
	// 0.125 ms
}

int getData(char port)
{
	unsigned char d[2];
	char t = 0;
	*d = 0x0000;

	PEOUT ^= 0x80;
	send_SPI(port); // Channel 3 - vRef
    t=100; while(--t);                      // Delay >5us
	read_SPI(&d);
	PEOUT = 0xFF;

	return (int)*d;
}

void test_SPI()
{
	char t = 0;
	unsigned char d[2];
	d[0] = 0x00;
	d[1] = 0x00;

	PEOUT ^= 0x80;
	send_SPI(0x1A); // Channel 3 - vRef
	waitForFullBuf();
    //t=0xFF; while(--t);                      // Delay >5us
	read_SPI(&d);
	writeData(d[0], d[1], d[0]);
	PEOUT = 0xFF;
}

/////////////////////////////////////////////////////////
/// Function: readData
///
/// <summary>Read data from the active accelerometer.</summary>
///
/// <returns>Nothing</returns>
/////////////////////////////////////////////////////////
void readData()
{
	int x = 0;
	int y = 0;
	int z = 0;
	float xg = 0.0;
	float yg = 0.0;
	float zg = 0.0;
	int vRef = 0;

	// 0.125 ms
	//t=2; while(--t); // 100ns

	vRef = getData(0x1B); // Channel 3 - vRef
	x = getData(0x18); // Channel 0
//	x = (x | MSB);
	xg = (x - vRef) * 0.0022;

	vRef = getData(0x1B); // Channel 3 - vRef
	y = getData(0x19); // Channel 1
//	y = (y | MSB);
	yg = (y - vRef) * 0.0022;

	vRef = getData(0x1B); // Channel 3 - vRef
	z = getData(0x1A); // Channel 2
//	z = (z | MSB);
	zg = (z - vRef) * 0.0022;

	// TODO : Transfer the captured data to the host pc.
	printf("%d %d %d %d\n", vRef, x, y, z);
//	writeData(x, y, z);
}

void initPortE()
{
	PEADDR = 0x02; // Port E Address reg selects Alternate
				   // Function sub-reg
	PECTL &= 0x00; // Port E Control reg set to no alternate
				   // function
	PEADDR = 0x01; // Port E Address reg selects
				   // Data Direction sub-reg
	PECTL &= 0x00; // Port E Control reg configures all bits
				   // as output
	PEADDR = 0x03; // Port E Address reg accessed
	 			   // Output Control sub-reg
	PECTL &= 0x00; // Port E Control reg sets Output Control
				   // to push-pull
}

/////////////////////////////////////////////////////////
/// Function: main
///
/// <summary></summary>
///
/// <returns>Nothing</returns>
/////////////////////////////////////////////////////////
int main ()
{
	int x = 0;
	int y = 0;
	int z = 0;

	initPortE();

	PEOUT = 0xFF; // for CS cuz it needs high for off - 0xFF;

	x = init_uart(_UART0,_DEFFREQ,_DEFBAUD); 
	if(x == 0)
		printf("Hello UART0\n"); // Write to _UART0

    Init_SPI();      // Master, Phase = 1, Clock polarity = 0
//	SPI_Init(0x8000, MASTER, TRUE, FALSE);      // Master, Phase = 1, Clock polarity = 0

	if(x == 0)
		printf("x, y, z\n");
	while(1)
	{
		PEOUT = 0xFF;
		test_SPI();
		PEOUT = 0xFF;
	    // TODO : cycle through each of the accelerometers and capture
		// the data
/*		int pin = 0x80; // Starts at PE7
		int i = 0;
		for(i = 0; i < g_NumAccel; i++)
		{
			PEOUT ^= pin;
			readData();
			PEOUT = 0xFF;
			pin = pin >> 1;
		}*/
	}
	return 0;
}


	// #define __DATA_DIRECTION          0x01
	// #define __ALTERNATE_FUNCTION      0x02
	// #define __OUTPUT_CONTROL          0x03
	// #define __HIGH_DRIVE_ENABLE       0x04
	// #define __SMR_ENABLE              0x05
	// PCADDR = __ALTERNATE_FUNCTION // Select Alternate Function
	// PCCTL = 0x36 // Enable SPI interface

	// PC3 SCK SPI Serial Clock
	// PC4 MOSI SPI Master Out Slave In
	// PC5 MISO SPI Master In Slave Out

/*
#ifdef EZ8_SPI
#define SPID    (*(unsigned char volatile far*)0xF60)              // Reset = 0xXX SPI Data
#define SPIDATA (*(unsigned char volatile far*)0xF60)              // Reset = 0xXX SPI Data
#define SPICTL  (*(unsigned char volatile far*)0xF61)              // Reset = 0x00 SPI Control
#define SPISTAT (*(unsigned char volatile far*)0xF62)              // Reset = 0x00 SPI Status
#define SPSTAT  (*(unsigned char volatile far*)0xF62)              // Reset = 0x00 SPI Status
#define SPIMODE (*(unsigned char volatile far*)0xF63)              // Reset = 0x00 SPI Mode

#if defined(_Z8F642) || defined(_Z8F08) || defined(_Z8FMC16)
#define SPIDST  (*(unsigned char volatile far*)0xF64)              // Reset = 0x00 SPI Diagnostic State
#endif

#define SPIBR   (*(unsigned int volatile far*)0xF66)               // Reset = 0xFFFF SPI Baud Rate
#define SPIBRH  (*(unsigned char volatile far*)0xF66)              // Reset = 0xFF SPI Baud Rate High
#define SPIBRL  (*(unsigned char volatile far*)0xF67)              // Reset = 0xFF SPI Baud Rate Low
#endif

#ifdef EZ8_ESPI
#define ESPIDATA  (*(unsigned char volatile far*)0xF60)              // Reset = 0xXX ESPI Data
#define ESPITDCR  (*(unsigned char volatile far*)0xF61)              // Reset = 0x00 ESPI Transmit Data Command
#define ESPICTL   (*(unsigned char volatile far*)0xF62)              // Reset = 0x00 ESPI Control
#define ESPIMODE  (*(unsigned char volatile far*)0xF63)              // Reset = 0x00 ESPI Mode
#define ESPISTAT  (*(unsigned char volatile far*)0xF64)              // Reset = 0x81 ESPI Status
#define ESPISTATE (*(unsigned char volatile far*)0xF65)              // Reset = 0x00 ESPI State
#define ESPIBR    (*(unsigned int volatile far*)0xF66)               // Reset = 0xFFFF ESPI Baud Rate
#define ESPIBRH   (*(unsigned char volatile far*)0xF66)              // Reset = 0xFF ESPI Baud Rate High
#define ESPIBRL   (*(unsigned char volatile far*)0xF67)              // Reset = 0xFF ESPI Baud Rate Low
#endif
*/