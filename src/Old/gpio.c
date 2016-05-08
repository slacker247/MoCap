/*************************************************
 *  Copyright (C) 1999-2008 by Zilog, Inc.
 *  All Rights Reserved
 *************************************************/

#include <eZ8.h>

/* Initializes LED ports - Port A */

void init_led_gpio(void)
{
  PAADDR 	= 0x01;     // PA Data Dir = output: updated
  PACTL 	= 0x00;    	// PA Out Ctrl = push-pull
  PAOUT		= 0X00;     // OUTPUT	           : updated
}

/* Initializes Test button port - Port C */

void init_test_button_gpio(void)
{
  PCADDR 	= 0x01;     // PC Data Dir = input: updated
  PCCTL 	= 0xFF;     // PC input Ctrl = Pin 0
}

/* Turns off ALL LED's */

void turn_off_led(void)
{
	PAOUT |= 0x07;  // Turn off all three leds
}


