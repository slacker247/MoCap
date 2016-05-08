// SPICTL Control Register Bits
#define PHASE           0x10    // Phase One, Data Read on falling edge of CLK
#define MMEN            0x02    // SPI Master Mode Enable
#define SPIEN           0x01    // SPI Enable
#define CLKPOL          0x08    // Clock Polarity

// SPIMODE Mode Register Bits
#define SSIO            0x02    // Slave Select I/O
#define SSV             0x01    // Slave Select Value

// SPISTAT Status Register Bits
#define TXST            0x02    // Transmit Status: 0=No Tx in progress, 1=Tx in progress

// BI2CBR Baud Rate Register Bits
#define Xtal            18432000ul      // Crystal frequency in Hz
#define SPIBaud         100ul           // Required SPI communications baud rate in kHz (100 kHz)
#define SPIBRG          Xtal/(2*SPIBaud*1000) // SPI baud rate register value

// General
#define SPIPort 		0x3C    // PC5=MISO, PC4=MOSI, PC3=SCK, PC2=SS
#define MSB             0x80    // To set or clear MSB in address

/*-------------------------------------------------------------------------------
**  Routine: Init_SPI
**  Parameters: void
**  Return: void
**  Purpose: Initalize the SPI controller to Master and set baud rate generator
--------------------------------------------------------------------------------*/
void Init_SPI(void);

void send_SPI(char data);
void read_SPI(unsigned char* data);
void waitForFullBuf();

/*-------------------------------------------------------------------------------
** Routine: ADEWrite
** Parameters: register address, data to write, no of bytes to write
** Return: void
** Purpose: Write data byte(s) to designated ADE7758 register
--------------------------------------------------------------------------------*/
void ADEWrite (unsigned char addr, unsigned char *data, unsigned char count);

/*-------------------------------------------------------------------------------
** Routine: ADERead
** Parameters: register address, data to read, no of bytes to read
** Return: None
** Purpose: Read data byte(s) from designated ADE7758 register
--------------------------------------------------------------------------------*/
void ADERead (unsigned char addr, unsigned char *data, unsigned char count);

/*-------------------------------------------------------------------------------
**  Routine: SPI_TDRE
**  Parameters: void
**  Return: void
**  Purpose: Poll waiting for the Transmit Data Register to empty
--------------------------------------------------------------------------------*/
static void SPI_TDRE (void);
static void SPI_RDRE (void);