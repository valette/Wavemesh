#ifndef rangedecod_h
#define rangedecod_h

/*
  rangedecoder.h     headerfile for range encoding

  (c) Michael Schindler
  1997, 1998, 1999, 2000
  http://www.compressconsult.com/
  michael@compressconsult.com

  C++ port by Sebastien Valette

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.  It may be that this
  program violates local patents in your country, however it is
  belived (NO WARRANTY!) to be patent-free here in Austria. Glen
  Langdon also confirmed my poinion that IBM UK did not protect that
  method.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.

  Range encoding is based on an article by G.N.N. Martin, submitted
  March 1979 and presented on the Video & Data Recording Conference,
  Southampton, July 24-27, 1979. If anyone can name the original
  copyright holder of that article please contact me; this might
  allow me to make that article available on the net for general
  public.

  Range coding is closely related to arithmetic coding, except that
  it does renormalisation in larger units than bits and is thus
  faster. An earlier version of this code was distributed as byte
  oriented arithmetic coding, but then I had no knowledge of Martin's
  paper from 1979.

  The input and output is done by the INBYTE and OUTBYTE macros
  defined in the .c file; change them as needed; the first parameter
  passed to them is a pointer to the rangecoder structure; extend that
  structure as needed (and don't forget to initialize the values in
  start_encoding resp. start_decoding). This distribution writes to
  stdout and reads from stdin.

  There are no global or static var's, so if the IO is thread save the
  whole rangecoder is - unless GLOBALRANGECODER in rangecod.h is defined.

  For error recovery the last 3 bytes written contain the total number
  of bytes written since starting the encoder. This can be used to
  locate the beginning of a block if you have only the end.
*/
#include <iostream>
#include <fstream>
#include <vtkObject.h>
#include "port.h"

extern char coderversion[];

typedef uint4 code_value;       /* Type of a rangecode value       */
                                /* must accomodate 32 bits          */
/* it is highly recommended that the total frequency count is less  */
/* than 1 << 19 to minimize rounding effects.                       */
/* the total frequency count MUST be less than 1<<23                */

typedef uint4 freq; 

class VTK_EXPORT rangedecoder
{
public:

		// open a file for coding
	void open_file(const char* file);

	// close the opened file
	void close_file();

/* Start the decoder                                         */
/* rc is the range coder to be used                          */
/* returns the char from start_encoding or EOF               */
int start_decoding();

/* Calculate culmulative frequency for next symbol. Does NO update!*/
/* rc is the range coder to be used                          */
/* tot_f is the total frequency                              */
/* or: totf is 1<<shift                                      */
/* returns the <= culmulative frequency                      */
freq decode_culfreq(freq tot_f );
freq decode_culshift(freq shift );

/* Update decoding state                                     */
/* rc is the range coder to be used                          */
/* sy_f is the interval length (frequency of the symbol)     */
/* lt_f is the lower end (frequency sum of < symbols)        */
/* tot_f is the total interval length (total frequency sum)  */
void decode_update(freq sy_f, freq lt_f, freq tot_f);
void decode_update_shift( freq f1, freq f2, freq f3) 
{this->decode_update((freq)(f1),(freq)(f2),(freq)1<<(f3));};

/* Decode a byte/short without modelling                     */
/* rc is the range coder to be used                          */
unsigned char decode_byte();
unsigned short decode_short();


/* Finish decoding                                           */
/* rc is the range coder to be used                          */
void done_decoding();

	rangedecoder(){};
	virtual ~rangedecoder(){};

protected:

	// this method is virtual so that class derivation allows other IO types
	virtual unsigned char inbyte();

private:
    uint4 low;           /* low end of interval */
    uint4 range;         /* length of interval */
    uint4 help;          /* bytes_to_follow resp. intermediate value */
    unsigned char buffer;/* buffer for input/output */

	std::fstream	InputFile;

void dec_normalize();

};
#endif
