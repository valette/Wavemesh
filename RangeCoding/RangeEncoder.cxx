/*
rangecod.c     range encoding

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
belived (NO WARRANTY!) to be patent-free. Glen Langdon also
confirmed my poinion that IBM UK did not protect that method.


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
paper from 1997.

The input and output is done by the inbyte and outbyte macros
defined in the .c file; change them as needed; the first parameter
passed to them is a pointer to the rangeencoder structure; extend that
structure as needed (and don't forget to initialize the values in
start_encoding resp. start_decoding). This distribution writes to
stdout and reads from stdin.

There are no global or static var's, so if the IO is thread save the
whole rangeencoder is - unless GLOBALRANGECODER in rangecod.h is defined.

For error recovery the last 3 bytes written contain the total number
of bytes written since starting the encoder. This can be used to
locate the beginning of a block if you have only the end. The header
size you pass to initrangecoder is included in that count.

There is a supplementary file called renorm95.c available at the
website (www.compressconsult.com/rangecoder/) that changes the range
coder to an arithmetic coder for speed comparisons.

define RENORM95 if you want the arithmetic coder type renormalisation.
Requires renorm95.c
Note that the old version does not write out the bytes since init.
you should not define GLOBALRANGECODER then. This Flag is provided
only for spped comparisons between both renormalizations, see my
data compression conference article 1998 for details.
*/
/* #define RENORM95 */

/*
define NOWARN if you do not expect more than 2^32 outstanding bytes 
since I recommend restarting the coder in intervals of less than    
2^23 symbols for error tolerance this is not expected
*/
#define NOWARN

/*
define EXTRAFAST for increased speed; you loose compression and
compatibility in exchange.
*/
/* #define EXTRAFAST */

#include "port.h"
#include "RangeEncoder.h"

/* SIZE OF RANGE ENCODING CODE VALUES. */

#define CODE_BITS 32
#define Top_value ((code_value)1 << (CODE_BITS-1))


#define SHIFT_BITS (CODE_BITS - 9)
#define EXTRA_BITS ((CODE_BITS-2) % 8 + 1)
#define Bottom_value (Top_value >> 8)

#ifdef NOWARN
char coderversion[]="rangecoder 1.3 NOWARN (c) 1997-2000 Michael Schindler";
#else   
char coderversion[]="rangecoder 1.3 (c) 1997-2000 Michael Schindler";
#endif   

void rangeencoder::outbyte(unsigned char a)
{
	unsigned char TempChar=a;
	this->OutputFile.write((char*)&TempChar, sizeof (char));
}

void rangeencoder::open_file(const char* file)
{ 
	this->OutputFile.open (file, std::ofstream::out | std::ofstream::trunc|std::ios::binary);
}
void rangeencoder::close_file()
{ 
	this->OutputFile.close();
};

/* rc is the range coder to be used                            */
/* c is written as first byte in the datastream                */
/* one could do without c, but then you have an additional if  */
/* per outputbyte.                                             */
void rangeencoder::start_encoding(char c, int initlength )
{  
	this->low = 0;                /* Full code range */
	this->range = Top_value;
	this->buffer = c;
	this->help = 0;               /* No bytes to follow */
	this->bytecount = initlength;
}

/* I do the normalization before I need a defined state instead of */
/* after messing it up. This simplifies starting and ending.       */
void rangeencoder::enc_normalize()
{   
	while(this->range <= Bottom_value)     /* do we need renormalisation?  */
	{ 
		if (this->low < (code_value)0xff<<SHIFT_BITS)  /* no carry possible --> output */
		{   outbyte(this->buffer);
		for(; this->help; this->help--)
			this->outbyte(0xff);
		this->buffer = (unsigned char)(this->low >> SHIFT_BITS);
		}
		else if (this->low & Top_value) /* carry now, no future carry */
		{   this->outbyte(this->buffer+1);
		for(; this->help; this->help--)
			this->outbyte(0);
		this->buffer = (unsigned char)(this->low >> SHIFT_BITS);
		} else                           /* passes on a potential carry */
#ifdef NOWARN
			this->help++;
#else
			if (this->bytestofollow++ == 0xffffffffL)
			{   		
				cerr<<"Too many bytes outstanding - File too large"<<endl;
				exit(1);
			}
#endif
			this->range <<= 8;
			this->low = (this->low<<8) & (Top_value-1);
			this->bytecount++;
	}
}

/* Encode a symbol using frequencies                         */
/* rc is the range coder to be used                          */
/* sy_f is the interval length (frequency of the symbol)     */
/* lt_f is the lower end (frequency sum of < symbols)        */
/* tot_f is the total interval length (total frequency sum)  */
/* or (faster): tot_f = (code_value)1<<shift                             */
void rangeencoder::encode_freq(freq sy_f, freq lt_f, freq tot_f )
{
	code_value r, tmp;
	this->enc_normalize();
	r = this->range / tot_f;
	tmp = r * lt_f;
	this->low += tmp;
#ifdef EXTRAFAST
	this->range = r * sy_f;
#else
	if (lt_f+sy_f < tot_f)
		this->range = r * sy_f;
	else
		this->range -= tmp;
#endif
}

void rangeencoder::encode_shift(freq sy_f, freq lt_f, freq shift )
{	
	code_value r, tmp;
	this->enc_normalize();
	r = this->range >> shift;
	tmp = r * lt_f;
	this->low += tmp;
#ifdef EXTRAFAST
	this->range = r * sy_f;
#else
	if ((lt_f+sy_f) >> shift)
		this->range -= tmp;
	else  
		this->range = r * sy_f;
#endif
}

/* Finish encoding                                           */
/* rc is the range coder to be used                          */
/* actually not that many bytes need to be output, but who   */
/* cares. I output them because decode will read them :)     */
/* the return value is the number of bytes written           */
uint4 rangeencoder::done_encoding()
{  
	uint tmp;
	this->enc_normalize();     /* now we have a normalized state */
	this->bytecount += 5;
	if ((this->low & (Bottom_value-1)) < ((this->bytecount&0xffffffL)>>1))
		tmp = this->low >> SHIFT_BITS;
	else
		tmp = (this->low >> SHIFT_BITS) + 1;
	if (tmp > 0xff) /* we have a carry */
	{   this->outbyte(this->buffer+1);
	for(; this->help; this->help--)
		this->outbyte(0);
	} else  /* no carry */
	{   this->outbyte(this->buffer);
	for(; this->help; this->help--)
		this->outbyte(0xff);
	}
	this->outbyte((tmp & 0xff));
	this->outbyte((this->bytecount>>16) & 0xff);
	this->outbyte((this->bytecount>>8) & 0xff);
	this->outbyte(this->bytecount & 0xff);
	return this->bytecount;
}
