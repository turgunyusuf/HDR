#include "vigra/impex.hxx"
//#include <vigra/filemapimage.hxx>
#include "vigra/colorconversions.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/functorexpression.hxx"
#include "HDRNumericTraits.h"
#include <map>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random.hpp>
#include <boost/assign/list_inserter.hpp>

#define MAX_PYRAMID_LEVELS 29  
#define DEFAULT_IMAGE_RESOLUTION 300.0f
#define NUMERIC_OPTION_DELIMITERS ";:/"

#define IMUL6(A) (A * SKIPSMImagePixelType(6))
#define IMUL5(A) (A * SKIPSMImagePixelType(5))
#define IMUL4(A) (A * SKIPSMImagePixelType(4))
#define IMUL11(A) (A * SKIPSMImagePixelType(11))
#define AMUL6(A) (A * SKIPSMAlphaPixelType(6))
using namespace vigra; 
using namespace std;

using vigra::Int8;
using vigra::Int16;
using vigra::Int32;
using vigra::Int64;
using vigra::UInt8;
using vigra::UInt16;
#ifndef WIN32
using vigra::UInt32;
#else
using vigra::UInt32;
#endif
using vigra::UInt64;
boost::mt19937 Twister;

// SKIPSM update routine used when visiting a pixel in the top two rows
// and the left two rows.
#define SKIPSM_EXPAND(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
	do {                                                            \
	current = SKIPSMImagePixelType(sa(sx));                     \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
	out01 = sc0a[srcx];                                         \
	out11 = sc0b[srcx];                                         \
	sc1a[srcx] = sc0a[srcx];                                    \
	sc1b[srcx] = sc0b[srcx];                                    \
	sc0a[srcx] = sr1 + IMUL6(sr0) + current;                    \
	sc0b[srcx] = (sr0 + current) * 4;                           \
	sr1 = sr0;                                                  \
	sr0 = current;                                              \
	out00 += sc0a[srcx];                                        \
	out10 += sc0b[srcx];                                        \
	out01 += sc0a[srcx];                                        \
	out11 += sc0b[srcx];                                        \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
	out11 /= SKIPSMImagePixelType(SCALE_OUT11);                 \
	da.set(cf(SKIPSMImagePixelType(da(dx)), out00), dx);        \
	++dx.x;                                                     \
	da.set(cf(SKIPSMImagePixelType(da(dx)), out10), dx);        \
	++dx.x;                                                     \
	da.set(cf(SKIPSMImagePixelType(da(dxx)), out01), dxx);      \
	++dxx.x;                                                    \
	da.set(cf(SKIPSMImagePixelType(da(dxx)), out11), dxx);      \
	++dxx.x;                                                    \
	} while (false)


// SKIPSM update routine used for the extra row under the main image body.
#define SKIPSM_EXPAND_ROW_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
	do {                                                            \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
	da.set(cf(da(dx), out00), dx);                              \
	++dx.x;                                                     \
	da.set(cf(da(dx), out10), dx);                              \
	++dx.x;                                                     \
	if (dst_h_even) {                                           \
	out01 = sc0a[srcx];                                     \
	out11 = sc0b[srcx];                                     \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);             \
	out11 /= SKIPSMImagePixelType(SCALE_OUT11);             \
	da.set(cf(da(dxx), out01), dxx);                        \
	++dxx.x;                                                \
	da.set(cf(da(dxx), out11), dxx);                        \
	++dxx.x;                                                \
	}                                                           \
	} while (false)


#define SKIPSM_EXPAND_COLUMN_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
	do {                                                            \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out01 = sc0a[srcx];                                         \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
	out11 = sc0b[srcx];                                         \
	sc1a[srcx] = sc0a[srcx];                                    \
	sc1b[srcx] = sc0b[srcx];                                    \
	sc0a[srcx] = sr1 + IMUL6(sr0);                              \
	sc0b[srcx] = sr0 * 4;                                       \
	out00 += sc0a[srcx];                                        \
	out01 += sc0a[srcx];                                        \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
	da.set(cf(da(dx), out00), dx);                              \
	da.set(cf(da(dxx), out01), dxx);                            \
	if (dst_w_even) {                                           \
	++dx.x;                                                 \
	++dxx.x;                                                \
	out10 += sc0b[srcx];                                    \
	out11 += sc0b[srcx];                                    \
	out10 /= SKIPSMImagePixelType(SCALE_OUT10);             \
	out11 /= SKIPSMImagePixelType(SCALE_OUT11);             \
	da.set(cf(da(dx), out10), dx);                          \
	da.set(cf(da(dxx), out11), dxx);                        \
	}                                                           \
	} while (false)


#define SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_EVEN(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
	do {                                                            \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out01 = sc0a[srcx];                                         \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
	out11 = sc0b[srcx];                                         \
	sc1a[srcx] = sc0a[srcx];                                    \
	sc1b[srcx] = sc0b[srcx];                                    \
	sc0a[srcx] = sr1 + IMUL6(sr0) + SKIPSMImagePixelType(sa(sy)); \
	sc0b[srcx] = (sr0 + SKIPSMImagePixelType(sa(sy))) * 4;      \
	out00 += sc0a[srcx];                                        \
	out01 += sc0a[srcx];                                        \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
	da.set(cf(da(dx), out00), dx);                              \
	da.set(cf(da(dxx), out01), dxx);                            \
	++dx.x;                                                     \
	++dxx.x;                                                    \
	out10 += sc0b[srcx];                                        \
	out11 += sc0b[srcx];                                        \
	out10 /= SKIPSMImagePixelType(SCALE_OUT10);                 \
	out11 /= SKIPSMImagePixelType(SCALE_OUT11);                 \
	da.set(cf(da(dx), out10), dx);                              \
	da.set(cf(da(dxx), out11), dxx);                            \
	} while (false)


#define SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_ODD(SCALE_OUT00, SCALE_OUT01) \
	do {                                                            \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out01 = sc0a[srcx];                                         \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                     \
	out11 = sc0b[srcx];                                         \
	sc1a[srcx] = sc0a[srcx];                                    \
	sc1b[srcx] = sc0b[srcx];                                    \
	sc0a[srcx] = sr1 + IMUL6(sr0) + IMUL4(SKIPSMImagePixelType(sa(sy))); \
	sc0b[srcx] = (sr0 + SKIPSMImagePixelType(sa(sy))) * 4;      \
	out00 += sc0a[srcx];                                        \
	out01 += sc0a[srcx];                                        \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);                 \
	da.set(cf(da(dx), out00), dx);                              \
	da.set(cf(da(dxx), out01), dxx);                            \
	} while (false)


#define SKIPSM_EXPAND_ROW_COLUMN_END(SCALE_OUT00, SCALE_OUT10, SCALE_OUT01, SCALE_OUT11) \
	do {                                                            \
	out00 = sc1a[srcx] + IMUL6(sc0a[srcx]);                     \
	out00 /= SKIPSMImagePixelType(SCALE_OUT00);                 \
	da.set(cf(da(dx), out00), dx);                              \
	if (dst_w_even) {                                           \
	out10 = sc1b[srcx] + IMUL6(sc0b[srcx]);                 \
	out10 /= SKIPSMImagePixelType(SCALE_OUT10);             \
	++dx.x;                                                 \
	da.set(cf(da(dx), out10), dx);                          \
	}                                                           \
	if (dst_h_even) {                                           \
	out01 = sc0a[srcx];                                     \
	out01 /= SKIPSMImagePixelType(SCALE_OUT01);             \
	da.set(cf(da(dxx), out01), dxx);                        \
	if (dst_w_even) {                                       \
	out11 = sc0b[srcx];                                 \
	out11 /= SKIPSMImagePixelType(SCALE_OUT11);         \
	++dxx.x;                                            \
	da.set(cf(da(dxx), out11), dxx);                    \
	}                                                       \
	}                                                           \
	} while (false)


typedef std::list<std::string> FileNameList;

//DEFINE_HDRNUMERICTRAITS(IMAGE_TYPE,vigra::UInt8,vigra::Int16)
//DEFINE_HDRNUMERICTRAITS(IMAGE_TYPE,vigra::UInt16,vigra::Int32)

struct ImageResolution 
{
	ImageResolution() : x(0.0f), y(0.0f) {}

	ImageResolution(float anXresolution, float aYresolution) :
	x(anXresolution), y(aYresolution) {}

	bool operator==(const ImageResolution& anOther) const 
	{
		return this->x == anOther.x && this->y == anOther.y;
	}

	bool operator!=(const ImageResolution& anOther) const
	{
		return !operator==(anOther);
	}

	float x;
	float y;
};

typedef enum BoundaryKind
{
	UnknownWrapAround,          // unknown kind
	OpenBoundaries,             // contractible
	HorizontalStrip,            // contractible along 2nd axis
	VerticalStrip,              // contractible along 1st axis
	DoubleStrip                 // non-contractible
} boundary_t;

struct AlternativePercentage
{
	double value;
	bool isPercentage;

	std::string str() const {
		std::ostringstream oss;
		oss << value;
		if (isPercentage) {oss << "%";}
		return oss.str();
	}
};

int ExactLevels = 0;
bool OneAtATime = true;
boundary_t WrapAround = OpenBoundaries;
bool GimpAssociatedAlphaHack_ = false;
bool UseCIECAM = false;
bool OutputSizeGiven = false;
int OutputWidthCmdLine = 0;
int OutputHeightCmdLine = 0;
int OutputOffsetXCmdLine = 0;
int OutputOffsetYCmdLine = 0;
std::string OutputCompression;
std::string OutputPixelType;
double WExposure = 1.0;         
double WContrast = 0.0;         
double WSaturation = 1.0;      
double WEntropy = 0.0;          
double WMu = 0.5;               
double WSigma = 0.2;           
bool WSaturationIsDefault = true;
int ContrastWindowSize = 5;
std::string GrayscaleProjector;
struct EdgeFilterConfiguration {double edgeScale, lceScale, lceFactor;} FilterConfig = {
	0.0,                       
	0.0,                      
	0.0                        
};
struct AlternativePercentage MinCurvature = {0.0, false}; 
int EntropyWindowSize = 3;      
struct AlternativePercentage EntropyLowerCutoff = {0.0, true}; 
struct AlternativePercentage EntropyUpperCutoff = {100.0, true}; 
bool UseHardMask = true;

typedef std::pair<double, double> range_t;

range_t rangeOfPixelType(const std::string& aPixelType)
{
	typedef std::map<std::string, range_t> Str2PairMapType;
	Str2PairMapType rangeMap;

	boost::assign::insert(rangeMap)
		("INT8", std::make_pair(vigra::NumericTraits<vigra::Int8>::min(),
		vigra::NumericTraits<vigra::Int8>::max()))
		("INT16", std::make_pair(vigra::NumericTraits<vigra::Int16>::min(),
		vigra::NumericTraits<vigra::Int16>::max()))
		("INT32", std::make_pair(vigra::NumericTraits<vigra::Int32>::min(),
		vigra::NumericTraits<vigra::Int32>::max()))

		("UINT8", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt8>::max()))
		("UINT16", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt16>::max()))
		("UINT32", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt32>::max()))

		("FLOAT", std::make_pair(0.0, 1.0))
		("DOUBLE", std::make_pair(0.0, 1.0));

	assert(!aPixelType.empty());
	Str2PairMapType::const_iterator r = rangeMap.find(aPixelType);
	if (r == rangeMap.end())
	{
		throw std::invalid_argument(std::string("unknown pixel type \"") + aPixelType + "\"");
	}
	else
	{
		return r->second;
	}
}
