#include <eZ8.h>
#include "SPI.h"

/*-------------------------------------------------------------------------------
**  Routine: Init_SPI
**  Parameters: void
**  Return: void
**  Purpose: Initalize the SPI controller to Master and set baud rate generator
--------------------------------------------------------------------------------*/
void Init_SPI(void)
{
        #define SPIPort 0x3C    // PC5=MISO, PC4=MOSI, PC3=SCK, PC2=SS

        // Set up SPI i/o pins 
        PCAF = SPIPort;                 // Alternate Function: Set PC5=MISO, PC4=MOSI, PC3=SCK, PC2=SS                                                  
        PCADDR = 0x00;                  // Advised to avoid inadvertent changes to port sub-registers

        // Set up SPI Registers
        SPICTL = PHASE | MMEN | SPIEN & ~CLKPOL;        // SCK Phase One, SPI Master, SPI enable, Clk idle low

        SPIMODE |= SSIO;                // Slave Select i/o: Set SS as o/p

        // Load the baud rate register to 100 kHz (ADE7758 can clock at 10 MHz)
        SPIBRH = (char)(SPIBRG >> 8);
        SPIBRL = (char)(SPIBRG & 0xFF);
}

/*-------------------------------------------------------------------------------
** Routine: ADEWrite
** Parameters: register address, data to write, no of bytes to write
** Return: void
** Purpose: Write data byte(s) to designated ADE7758 register
--------------------------------------------------------------------------------*/
void ADEWrite (unsigned char addr, unsigned char *data, unsigned char count)
{
        unsigned char i;
        unsigned char t;

        DI();                                           // Disable interrupts
        SPIMODE &= ~SSV;                        // Slave Select: Clear SS low to select ADE7758
        SPIDATA = (addr | MSB);         // Write address to SPI Data Register, MSB =1
        SPI_TDRE();                                     // Wait for Transmit Buffer Empty
       
        t=2; while(--t);                        // Delay 100ns
       
        for (i=0; i<count; i++)
        {
                SPIDATA = *data--;              // Write data byte(s) to SPI Data Register
                SPI_TDRE();                             // Wait for Transmit Buffer Empty
        }

        SPIMODE |= SSV;                         // Slave Select: Set SS high to deselect ADE7758
        EI();                                           // Re-establish interrupts
}

/*-------------------------------------------------------------------------------
** Routine: ADERead
** Parameters: register address, data to read, no of bytes to read
** Return: None
** Purpose: Read data byte(s) from designated ADE7758 register
--------------------------------------------------------------------------------*/
void ADERead (unsigned char addr, unsigned char *data, unsigned char count)
{
        unsigned char i;
        unsigned char t;

        DI();                                           // Disable interrupts
        SPIMODE &= ~SSV;                        // Slave Select: Clear SS low to select ADE7758
        SPIDATA = (addr & ~MSB);        // Write Address to SPI Data Reg, MSB =0.
        SPI_TDRE();                                     // Wait for Transmit Buffer Empty
       
        t=100; while(--t);                      // Delay >5us

        for (i=0; i<count; i++)
        {
                SPIDATA = 0xF0;                 // Dummy write forces data onto DOUT line (into MISO)
                SPI_TDRE();                             // Wait for Transmit Buffer Empty
                *data-- = SPIDATA;              // Read data byte(s) from SPI Data Register
        }
       
        SPIMODE |= SSV;                         // Slave Select: Set SS high to deselect ADE7758
        EI();                                           // Re-establish interrupts ??
}

/*-------------------------------------------------------------------------------
**  Routine: SPI_TDRE
**  Parameters: void
**  Return: void
**  Purpose: Poll waiting for the Transmit Data Register to empty
--------------------------------------------------------------------------------*/
static void SPI_TDRE (void)
{
        while ((SPISTAT & TXST) == 0x00);       //Wait for Transmit Buffer Full
        while ((SPISTAT & TXST) == TXST);       //Wait for Transmit Buffer Empty
}

/* --- General Include Files --- */
#include <eZ8.h>
#include <stdio.h>
//#include "MAIN.h"                             // Main program

#include "SPI.h"                                // Serial peripherial interface

/* --- Global Variables --- */
unsigned int volatile current_ms = 0;

//-------------------------------- Main Loop ---------------------------------
void main()
{
        unsigned char i;
        unsigned char Rx=0;
       
        /* --- Initialising routines (don't change order) --- */
        EI();                           // Master enable interrupts
        Init_T0();                      // Initialise ms timer
        Init_SPI();                     // Initialise the serial interface to the ADE7758

        /* --- Write splash message to UART0 --- */
        printf("Test SPI to ADE7758\n");

 
        /* --- Test ADE7758 --- */
        i = 0x44;                                               // SWRST bit + DISCF bit 
        ADEWrite (OPMODE, &i, 1);
        Delayms (1);                                    // Need to wait about 500us after reset
//      ADERead (RSTATUS,(unsigned char *) &j, 2);              // Read/clear STATUS after SWRST 
 
        i = 0x02;
        ADEWrite (WAVMODE, &i, 1);

        ADERead (WAVMODE, &Rx, 1);              // Get back 0x03 so wrong!

        ADERead (VERS, &Rx, 1);                 // VERSION of Die [expect 0x81 but get 0x03] 
        printf("DIE = %#2x \n", Rx);
}


