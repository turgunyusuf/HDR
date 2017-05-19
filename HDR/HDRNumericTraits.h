#ifndef _HDR_NUMERIC_TRAITS_H_
#define _HDR_NUMERIC_TRAITS_H_

#include "vigra\filemapimage.hxx"
#include "vigra\basicimage.hxx"

using namespace vigra;
#define IMAGE_TYPE FileMapImage
struct Error_hdrNumericTraits_not_specialized_for_this_case { };

template<class A>

struct hdrNumericTraits 
{
	//定义输入原始图像的像素类型
	typedef Error_hdrNumericTraits_not_specialized_for_this_case PixelType;


	//定义图像模板类型
    typedef Error_hdrNumericTraits_not_specialized_for_this_case MaskImage;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case AlphaImage;
	
	enum { ImagePyramidIntegerBits = 0 };
	enum { ImagePyramidFractionBits = 0 };

	enum { MaskPyramidIntegerBits = 0 };
	enum { MaskPyramidFractionBits = 0 };

	//定义金字塔图像类型
	typedef  Error_hdrNumericTraits_not_specialized_for_this_case PyramidImage;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case ImagePyramidPixelComponentType;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case ImagePyramidPixelType;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case ImagePyramidType;
	// Pixel type used by SKIPSM algorithm for intermediate image pixel calculations
	typedef Error_hdrNumericTraits_not_specialized_for_this_case SKIPSMImagePixelType;

	// Pixel type used by SKIPSM algorithm for intermediate alpha pixel calculations
	typedef Error_hdrNumericTraits_not_specialized_for_this_case SKIPSMAlphaPixelType;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case SKIPSMMaskPixelType;
	typedef Error_hdrNumericTraits_not_specialized_for_this_case MaskPyramidType;
};

#define DEFINE_HDRNUMERICTRAITS(IMAGE_TYPE, IMAGECOMPONENT,PYRAMIDCOMPONENT, PYRAMIDINTEGER,PYRAMIDFRACTION,MASKPYRAMIDINTEGER,MASKPYRAMIDFRACTION,MASKPYRAMID)\
	\
  template<>\
struct hdrNumericTraits<IMAGECOMPONENT>{\
	typedef IMAGECOMPONENT PixelType; \
	typedef IMAGE_TYPE<IMAGECOMPONENT> ImageType; \
	typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType; \
	typedef PYRAMIDCOMPONENT ImagePyramidPixelType; \
	typedef IMAGE_TYPE<PYRAMIDCOMPONENT> ImagePyramidType; \
	enum {ImagePyramidIntegerBits = PYRAMIDINTEGER}; \
	enum {ImagePyramidFractionBits = PYRAMIDFRACTION}; \
	enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER}; \
	enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION}; \
	typedef MASKPYRAMID MaskPyramidPixelType; \
	typedef IMAGE_TYPE<MASKPYRAMID> MaskPyramidType; \
};\
	template<> \
struct hdrNumericTraits<RGBValue<IMAGECOMPONENT,0,1,2> > { \
	typedef RGBValue<IMAGECOMPONENT,0,1,2> PixelType; \
	typedef IMAGE_TYPE<RGBValue<IMAGECOMPONENT,0,1,2> > ImageType; \
	typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType; \
	typedef RGBValue<PYRAMIDCOMPONENT,0,1,2> ImagePyramidPixelType; \
	typedef IMAGE_TYPE<RGBValue<PYRAMIDCOMPONENT,0,1,2> > ImagePyramidType; \
	enum {ImagePyramidIntegerBits = PYRAMIDINTEGER}; \
	enum {ImagePyramidFractionBits = PYRAMIDFRACTION}; \
	enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER}; \
	enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION}; \
	typedef MASKPYRAMID MaskPyramidPixelType; \
	typedef IMAGE_TYPE<MASKPYRAMID> MaskPyramidType; \
};

  DEFINE_HDRNUMERICTRAITS(IMAGE_TYPE, UInt8,Int16,9,7,9,7,Int16)
  DEFINE_HDRNUMERICTRAITS(IMAGE_TYPE, UInt16,Int32,17,7,9,15,Int32)



#define ALPHA_TRAITS(T1,S) \
	template<> \
   struct AlphaTraits<T1> \
   { \
   static T1 max() \
   { \
   return S; \
   } \
   static T1 zero() \
   { \
   return 0; \
   } \
   }; \
   template<> \
   struct AlphaTraits<vigra::RGBValue<T1> > \
   { \
   static T1 max() \
   { \
   return S; \
   } \
   static T1 zero() \
   { \
   return 0; \
   } \
   };

	template <class T1>
   struct AlphaTraits;

   ALPHA_TRAITS(unsigned char, UCHAR_MAX);
   ALPHA_TRAITS(signed char, SCHAR_MAX);
   ALPHA_TRAITS(unsigned short, USHRT_MAX);
   ALPHA_TRAITS(signed short, SHRT_MAX);
   ALPHA_TRAITS(unsigned int, UINT_MAX);
   ALPHA_TRAITS(signed int, INT_MAX);
   ALPHA_TRAITS(float, 1.0);
   ALPHA_TRAITS(double, 1.0);

#endif