#ifndef rangedecod2_h
#define rangedecod2_h

#include "RangeDecoder.h"

class VTK_EXPORT RangeDecoderHook
{
public:
	virtual void ByteRead(){};
	
	RangeDecoderHook(){};
	~RangeDecoderHook(){};
};

class rangedecoder2 : public rangedecoder
{
public:

	void SetHook(RangeDecoderHook *Hook)
	{
		this->Hook=Hook;
	}
	
	rangedecoder2()
	{
		this->Hook=0;
	};
	
	~rangedecoder2(){};

protected:

	virtual unsigned char inbyte()
	{
		if (this->Hook)
			this->Hook->ByteRead();

		return this->rangedecoder::inbyte();
	};
	
	RangeDecoderHook *Hook;
};
#endif
