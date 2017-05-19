#include "common.h"
//#include <vigra/functorexpression.hxx>
//#include "openmp.h"
static inline double square(double x)
{
	return x * x;
}


template <typename SKIPSMImagePixelType, typename PyramidImageType>
void
	collapsePyramid(bool wraparound, vector<PyramidImageType*> *p)
{

	for (int l = (p->size()-2); l >= 0; l--)
	{ 
		expand<SKIPSMImagePixelType>(true, wraparound,
			srcImageRange(*((*p)[l + 1])),
			destImageRange(*((*p)[l])));
	}

}


template <typename MaskPixelType>
class ImageMaskMultiplyFunctor {
public:
	ImageMaskMultiplyFunctor(MaskPixelType d) :
	  divisor(NumericTraits<MaskPixelType>::toRealPromote(d)) {}

	  template <typename ImagePixelType>
	  ImagePixelType operator()(const ImagePixelType& iP, const MaskPixelType& maskP ) const 
	  {
		  typedef typename NumericTraits<ImagePixelType>::RealPromote RealImagePixelType;

		  const double maskCoeff = NumericTraits<MaskPixelType>::toRealPromote(maskP) / divisor;
		  const RealImagePixelType riP = NumericTraits<ImagePixelType>::toRealPromote(iP);
		  const RealImagePixelType blendP = riP * maskCoeff;
		  return NumericTraits<ImagePixelType>::fromRealPromote(blendP);
	  }

protected:
	double divisor;
};


template <typename SKIPSMImagePixelType,
	typename SrcImageIterator, typename SrcAccessor,
	typename DestImageIterator, typename DestAccessor,
	typename CombineFunctor>
	void
	expand(bool add, bool wraparound,
	SrcImageIterator src_upperleft,
	SrcImageIterator src_lowerright,
	SrcAccessor sa,
	DestImageIterator dest_upperleft,
	DestImageIterator dest_lowerright,
	DestAccessor da,
	CombineFunctor cf)
{
	int src_w = src_lowerright.x - src_upperleft.x;
	int src_h = src_lowerright.y - src_upperleft.y;
	int dst_w = dest_lowerright.x - dest_upperleft.x;
	int dst_h = dest_lowerright.y - dest_upperleft.y;

	const bool dst_w_even = (dst_w & 1) == 0;
	const bool dst_h_even = (dst_h & 1) == 0;

	SKIPSMImagePixelType current;
	SKIPSMImagePixelType out00, out10, out01, out11;
	SKIPSMImagePixelType sr0, sr1;
	SKIPSMImagePixelType* sc0a = new SKIPSMImagePixelType[src_w + 1];
	SKIPSMImagePixelType* sc0b = new SKIPSMImagePixelType[src_w + 1];
	SKIPSMImagePixelType* sc1a = new SKIPSMImagePixelType[src_w + 1];
	SKIPSMImagePixelType* sc1b = new SKIPSMImagePixelType[src_w + 1];

	const SKIPSMImagePixelType SKIPSMImageZero(NumericTraits<SKIPSMImagePixelType>::zero());

	DestImageIterator dy = dest_upperleft;
	DestImageIterator dyy = dest_upperleft;
	DestImageIterator dx = dy;
	DestImageIterator dxx = dyy;
	SrcImageIterator sy = src_upperleft;
	SrcImageIterator sx = sy;

	int srcy = 0;
	int srcx = 0;
	{
		srcx = 0;
		sx = sy;
		sr0 = SKIPSMImagePixelType(sa(sx));
		if (wraparound) {
			sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 1, 0)));
			if (!dst_w_even) {
				sr1 = IMUL4(sr1);
			}
		} else {
			sr1 = SKIPSMImageZero;
		}

		srcx = 1;
		++sx.x;

		for (; srcx < src_w; ++srcx, ++sx.x) {
			current = SKIPSMImagePixelType(sa(sx));
			sc0a[srcx] = sr1 + IMUL6(sr0) + current;
			sc0b[srcx] = (sr0 + current) * 4;
			sc1a[srcx] = SKIPSMImageZero;
			sc1b[srcx] = SKIPSMImageZero;
			sr1 = sr0;
			sr0 = current;
		}

		if (wraparound) {
			current = SKIPSMImagePixelType(sa(sy));
			if (dst_w_even) {
				sc0a[srcx] = sr1 + IMUL6(sr0) + current;
				sc0b[srcx] = (sr0 + current) * 4;
			} else {
				sc0a[srcx] = sr1 + IMUL6(sr0) + IMUL4(current);
			}
		} else {
			sc0a[srcx] = sr1 + IMUL6(sr0);
			sc0b[srcx] = sr0 * 4;
		}
		sc1a[srcx] = SKIPSMImageZero;
		sc1b[srcx] = SKIPSMImageZero;
	}

	++dyy.y;
	srcy = 1;
	++sy.y;

	if (src_h > 1) {
		srcx = 0;
		sx = sy;
		sr0 = SKIPSMImagePixelType(sa(sx));
		if (wraparound) {
			sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 1, 0)));
			if (!dst_w_even) {
				sr1 = IMUL4(sr1);
			}
		} else {
			sr1 = SKIPSMImageZero;
		}

		srcx = 1;
		++sx.x;
		dx = dy;
		dxx = dyy;

		if (src_w > 1) {
			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND(56, 56, 16, 16);
				} else {
					SKIPSM_EXPAND(77, 56, 22, 16);
				}
			} else {
				SKIPSM_EXPAND(49, 56, 14, 16);
			}

			for (srcx = 2, ++sx.x; srcx < src_w; ++srcx, ++sx.x) {
				SKIPSM_EXPAND(56, 56, 16, 16);
			}

			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_EVEN(56, 56, 16, 16);
				} else {
					SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_ODD(77, 22);
				}
			} else {
				SKIPSM_EXPAND_COLUMN_END(49, 28, 14, 8);
			}
		} else {
			SKIPSM_EXPAND_COLUMN_END(42, 28, 12, 8);
		}
	} else {
		srcx = 0;
		sr0 = SKIPSMImageZero;
		sr1 = SKIPSMImageZero;

		dx = dy;
		dxx = dyy;

		if (src_w > 1) {
			srcx = 1;
			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_ROW_END(48, 48, 8, 8);
				} else {
					SKIPSM_EXPAND_ROW_END(66, 48, 11, 8);
				}
			} else {
				SKIPSM_EXPAND_ROW_END(42, 48, 7, 8);
			}

			for (srcx = 2; srcx < src_w; ++srcx) {
				SKIPSM_EXPAND_ROW_END(48, 48, 8, 8);
			}

			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_ROW_COLUMN_END(48, 48, 8, 8);
				} else {
					SKIPSM_EXPAND_ROW_COLUMN_END(66, 48, 11, 8);
				}
			} else {
				SKIPSM_EXPAND_ROW_COLUMN_END(42, 24, 7, 4);
			}
		} else {
			SKIPSM_EXPAND_ROW_COLUMN_END(36, 24, 6, 4);
		}

		delete [] sc0a;
		delete [] sc0b;
		delete [] sc1a;
		delete [] sc1b;

		return;
	}

	dy.y += 2;
	dyy.y += 2;
	srcy = 2;
	++sy.y;

	for (srcy = 2, sx = sy; srcy < src_h; ++srcy, ++sy.y, dy.y += 2, dyy.y += 2) {
		srcx = 0;
		sx = sy;
		sr0 = SKIPSMImagePixelType(sa(sx));
		if (wraparound) {
			sr1 = SKIPSMImagePixelType(sa(sy, Diff2D(src_w-1,0)));
			if (!dst_w_even) {
				sr1 = IMUL4(sr1);
			}
		} else {
			sr1 = SKIPSMImageZero;
		}

		srcx = 1;
		++sx.x;
		dx = dy;
		dxx = dyy;

		if (src_w > 1) {
			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND(64, 64, 16, 16);
				} else {
					SKIPSM_EXPAND(88, 64, 22, 16);
				}
			} else {
				SKIPSM_EXPAND(56, 64, 14, 16);
			}

			for (srcx = 2, ++sx.x; srcx < src_w; ++srcx, ++sx.x)
			{
				SKIPSM_EXPAND(64, 64, 16, 16);
			}

			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_EVEN(64, 64, 16, 16);
				} else {
					SKIPSM_EXPAND_COLUMN_END_WRAPAROUND_ODD(88, 22);
				}
			} else {
				SKIPSM_EXPAND_COLUMN_END(56, 32, 14, 8);
			}
		}
		else
		{
			SKIPSM_EXPAND_COLUMN_END(48, 32, 12, 8);
		}
	}

	{
		srcx = 0;
		sr0 = SKIPSMImageZero;
		sr1 = SKIPSMImageZero;

		dx = dy;
		dxx = dyy;

		if (src_w > 1)
		{
			srcx = 1;
			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_ROW_END(56, 56, 8, 8);
				} else {
					SKIPSM_EXPAND_ROW_END(77, 56, 11, 8);
				}
			} else {
				SKIPSM_EXPAND_ROW_END(49, 56, 7, 8);
			}

			for (srcx = 2; srcx < src_w; ++srcx) {
				SKIPSM_EXPAND_ROW_END(56, 56, 8, 8);
			}

			if (wraparound) {
				if (dst_w_even) {
					SKIPSM_EXPAND_ROW_COLUMN_END(56, 56, 8, 8);
				} else {
					SKIPSM_EXPAND_ROW_COLUMN_END(77, 56, 11, 8);
				}
			} else {
				SKIPSM_EXPAND_ROW_COLUMN_END(49, 28, 7, 4);
			}
		} 
		else 
		{
			SKIPSM_EXPAND_ROW_COLUMN_END(42, 28, 6, 4);
		}
	}

	delete [] sc0a;
	delete [] sc0b;
	delete [] sc1a;
	delete [] sc1b;
}

template<typename T1, typename T2, typename T3>
struct FromPromotePlusFunctorWrapper :
	public std::binary_function<T1, T2, T3>
{
	inline T3 operator()(const T1& a, const T2& b) const 
	{
		return NumericTraits<T3>::fromPromote(a + b);
	}
};

template <typename SKIPSMImagePixelType,
	typename SrcImageIterator, typename SrcAccessor,
	typename DestImageIterator, typename DestAccessor>
	inline void
	expand(bool add, bool wraparound,
	triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
	triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
{
	typedef typename DestAccessor::value_type DestPixelType;


	if (add)
	{
		expand<SKIPSMImagePixelType>(add, wraparound,
			src.first, src.second, src.third,
			dest.first, dest.second, dest.third,
			FromPromotePlusFunctorWrapper<DestPixelType, SKIPSMImagePixelType, DestPixelType>());
	}
	else
	{
		expand<SKIPSMImagePixelType>(add, wraparound,
			src.first, src.second, src.third,
			dest.first, dest.second, dest.third,
			std::minus<SKIPSMImagePixelType>());
	}

}



template <typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType,
	typename SrcImageIterator, typename SrcAccessor,
	typename AlphaIterator, typename AlphaAccessor,
	typename DestImageIterator, typename DestAccessor,
	typename DestAlphaIterator, typename DestAlphaAccessor>
	inline void
	reduce(bool wraparound,
	SrcImageIterator src_upperleft,
	SrcImageIterator src_lowerright,
	SrcAccessor sa,
	AlphaIterator alpha_upperleft,
	AlphaAccessor aa,
	DestImageIterator dest_upperleft,
	DestImageIterator dest_lowerright,
	DestAccessor da,
	DestAlphaIterator dest_alpha_upperleft,
	DestAlphaIterator dest_alpha_lowerright,
	DestAlphaAccessor daa)
{
	typedef typename DestAccessor::value_type DestPixelType;
	typedef typename DestAlphaAccessor::value_type DestAlphaPixelType;

	int src_w = src_lowerright.x - src_upperleft.x;
	int src_h = src_lowerright.y - src_upperleft.y;
	int dst_w = dest_lowerright.x - dest_upperleft.x;
	//int dst_h = dest_lowerright.y - dest_upperleft.y;

	vigra_precondition(src_w > 1 && src_h > 1,
		"src image too small in reduce");

	// State variables for source image pixel values
	SKIPSMImagePixelType isr0, isr1, isrp;
	SKIPSMImagePixelType* isc0 = new SKIPSMImagePixelType[dst_w + 1];
	SKIPSMImagePixelType* isc1 = new SKIPSMImagePixelType[dst_w + 1];
	SKIPSMImagePixelType* iscp = new SKIPSMImagePixelType[dst_w + 1];

	// State variables for source mask pixel values
	SKIPSMAlphaPixelType asr0, asr1, asrp;
	SKIPSMAlphaPixelType* asc0 = new SKIPSMAlphaPixelType[dst_w + 1];
	SKIPSMAlphaPixelType* asc1 = new SKIPSMAlphaPixelType[dst_w + 1];
	SKIPSMAlphaPixelType* ascp = new SKIPSMAlphaPixelType[dst_w + 1];

	// Convenient constants
	const SKIPSMImagePixelType SKIPSMImageZero(NumericTraits<SKIPSMImagePixelType>::zero());
	const SKIPSMAlphaPixelType SKIPSMAlphaZero(NumericTraits<SKIPSMAlphaPixelType>::zero());
	const SKIPSMAlphaPixelType SKIPSMAlphaOne(NumericTraits<SKIPSMAlphaPixelType>::one());
	const DestPixelType DestImageZero(NumericTraits<DestPixelType>::zero());
	const DestAlphaPixelType DestAlphaZero(NumericTraits<DestAlphaPixelType>::zero());
	const DestAlphaPixelType DestAlphaMax(NumericTraits<DestAlphaPixelType>::max());

	DestImageIterator dy = dest_upperleft;
	DestImageIterator dx = dy;
	SrcImageIterator sy = src_upperleft;
	SrcImageIterator sx = sy;
	AlphaIterator ay = alpha_upperleft;
	AlphaIterator ax = ay;
	DestAlphaIterator day = dest_alpha_upperleft;
	DestAlphaIterator dax = day;

	bool evenY = true;
	bool evenX = true;
	int srcy = 0;
	int srcx = 0;
	//int dsty = 0;
	int dstx = 0;

	// First row
	{
		if (wraparound) {
			asr0 = aa(ay, Diff2D(src_w - 2, 0)) ? SKIPSMAlphaOne : SKIPSMAlphaZero;
			asr1 = SKIPSMAlphaZero;
			asrp = aa(ay, Diff2D(src_w - 1, 0)) ? (SKIPSMAlphaOne * 4) : SKIPSMAlphaZero;
			isr0 = aa(ay, Diff2D(src_w - 2, 0)) ? SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 2, 0))) : SKIPSMImageZero;
			isr1 = SKIPSMImageZero;
			isrp = aa(ay, Diff2D(src_w - 1, 0)) ? (SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 1, 0))) * 4) : SKIPSMImageZero;
		} else {
			asr0 = SKIPSMAlphaZero;
			asr1 = SKIPSMAlphaZero;
			asrp = SKIPSMAlphaZero;
			isr0 = SKIPSMImageZero;
			isr1 = SKIPSMImageZero;
			isrp = SKIPSMImageZero;
		}

		for (sx = sy, ax = ay, evenX = true, srcx = 0, dstx = 0;
			srcx < src_w;
			++srcx, ++sx.x, ++ax.x) {
				SKIPSMAlphaPixelType mcurrent(aa(ax) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
				SKIPSMImagePixelType icurrent(aa(ax) ? SKIPSMImagePixelType(sa(sx)) : SKIPSMImageZero);
				if (evenX) {
					asc1[dstx] = SKIPSMAlphaZero;
					asc0[dstx] = asr1 + AMUL6(asr0) + asrp + mcurrent;
					asr1 = asr0 + asrp;
					asr0 = mcurrent;
					isc1[dstx] = SKIPSMImageZero;
					isc0[dstx] = isr1 + IMUL6(isr0) + isrp + icurrent;
					isr1 = isr0 + isrp;
					isr0 = icurrent;
				} else {
					asrp = mcurrent * 4;
					isrp = icurrent * 4;
					++dstx;
				}
				evenX = !evenX;
		}

		// Last entries in first row
		if (!evenX) {
			// previous srcx was even
			++dstx;
			if (wraparound) {
				asc1[dstx] = SKIPSMAlphaZero;
				asc0[dstx] = asr1 + AMUL6(asr0) + (aa(ay) ? (SKIPSMAlphaOne * 4) : SKIPSMAlphaZero) + (aa(ay, Diff2D(1,0)) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
				isc1[dstx] = SKIPSMImageZero;
				isc0[dstx] = isr1 + IMUL6(isr0) + (aa(ay) ? (SKIPSMImagePixelType(sa(sy)) * 4) : SKIPSMImageZero)
					+ (aa(ay, Diff2D(1, 0)) ? SKIPSMImagePixelType(sa(sy, Diff2D(1, 0))) : SKIPSMImageZero);
			} else {
				asc1[dstx] = SKIPSMAlphaZero;
				asc0[dstx] = asr1 + AMUL6(asr0);
				isc1[dstx] = SKIPSMImageZero;
				isc0[dstx] = isr1 + IMUL6(isr0);
			}
		} else {
			// previous srcx was odd
			if (wraparound) {
				asc1[dstx] = SKIPSMAlphaZero;
				asc0[dstx] = asr1 + AMUL6(asr0) + asrp + (aa(ay) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
				isc1[dstx] = SKIPSMImageZero;
				isc0[dstx] = isr1 + IMUL6(isr0) + isrp + (aa(ay) ? SKIPSMImagePixelType(sa(sy)) : SKIPSMImageZero);
			} else {
				asc1[dstx] = SKIPSMAlphaZero;
				asc0[dstx] = asr1 + AMUL6(asr0) + asrp;
				isc1[dstx] = SKIPSMImageZero;
				isc0[dstx] = isr1 + IMUL6(isr0) + isrp;
			}
		}
	}
	++sy.y;
	++ay.y;

	// Main Rows
	{
		for (evenY = false, srcy = 1; srcy < src_h; ++srcy, ++sy.y, ++ay.y) {
			if (wraparound) {
				asr0 = aa(ay, Diff2D(src_w - 2, 0)) ? SKIPSMAlphaOne : SKIPSMAlphaZero;
				asr1 = SKIPSMAlphaZero;
				asrp = aa(ay, Diff2D(src_w - 1, 0)) ? (SKIPSMAlphaOne * 4) : SKIPSMAlphaZero;
				isr0 = aa(ay, Diff2D(src_w - 2, 0)) ? SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 2,0))) : SKIPSMImageZero;
				isr1 = SKIPSMImageZero;
				isrp = aa(ay, Diff2D(src_w - 1, 0)) ? (SKIPSMImagePixelType(sa(sy, Diff2D(src_w - 1,0))) * 4) : SKIPSMImageZero;
			} else {
				asr0 = SKIPSMAlphaZero;
				asr1 = SKIPSMAlphaZero;
				asrp = SKIPSMAlphaZero;
				isr0 = SKIPSMImageZero;
				isr1 = SKIPSMImageZero;
				isrp = SKIPSMImageZero;
			}

			if (evenY) {
				// Even-numbered row

				// First entry in row
				sx = sy;
				ax = ay;
				if (wraparound) {
					asr1 = asr0 + asrp;
					isr1 = isr0 + isrp;
				}
				asr0 = aa(ax) ? SKIPSMAlphaOne : SKIPSMAlphaZero;
				isr0 = aa(ax) ? SKIPSMImagePixelType(sa(sx)) : SKIPSMImageZero;
				// isc*[0] are never used
				++sx.x;
				++ax.x;
				dx = dy;
				dax = day;

				// Main entries in row
				for (evenX = false, srcx = 1, dstx = 0; srcx < src_w; ++srcx, ++sx.x, ++ax.x) {
					SKIPSMAlphaPixelType mcurrent(aa(ax) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
					SKIPSMImagePixelType icurrent(aa(ax) ? SKIPSMImagePixelType(sa(sx)) : SKIPSMImageZero);
					if (evenX) {
						SKIPSMAlphaPixelType ap = asc1[dstx] + AMUL6(asc0[dstx]) + ascp[dstx];
						asc1[dstx] = asc0[dstx] + ascp[dstx];
						asc0[dstx] = asr1 + AMUL6(asr0) + asrp + mcurrent;
						asr1 = asr0 + asrp;
						asr0 = mcurrent;
						ap += asc0[dstx];

						SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
						isc1[dstx] = isc0[dstx] + iscp[dstx];
						isc0[dstx] = isr1 + IMUL6(isr0) + isrp + icurrent;
						isr1 = isr0 + isrp;
						isr0 = icurrent;
						if (ap) {
							ip += isc0[dstx];
							ip /= SKIPSMImagePixelType(ap);
							da.set(DestPixelType(ip), dx);
							daa.set(DestAlphaMax, dax);
						} else {
							da.set(DestImageZero, dx);
							daa.set(DestAlphaZero, dax);
						}

						++dx.x;
						++dax.x;
					} else {
						asrp = mcurrent * 4;
						isrp = icurrent * 4;
						++dstx;
					}
					evenX = !evenX;
				}

				// Last entries in row
				if (!evenX) {
					// previous srcx was even
					++dstx;

					SKIPSMAlphaPixelType ap = asc1[dstx] + AMUL6(asc0[dstx]) + ascp[dstx];
					asc1[dstx] = asc0[dstx] + ascp[dstx];
					if (wraparound) {
						asc0[dstx] = asr1 + AMUL6(asr0) + (aa(ay) ? (SKIPSMAlphaOne * 4) : SKIPSMAlphaZero) + (aa(ay, Diff2D(1,0)) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
					} else {
						asc0[dstx] = asr1 + AMUL6(asr0);
					}
					ap += asc0[dstx];

					SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
					isc1[dstx] = isc0[dstx] + iscp[dstx];
					if (wraparound) {
						isc0[dstx] = isr1 + IMUL6(isr0) + (aa(ay) ? (SKIPSMImagePixelType(sa(sy)) * 4) : SKIPSMImageZero)
							+ (aa(ay, Diff2D(1, 0)) ? SKIPSMImagePixelType(sa(sy, Diff2D(1, 0))) : SKIPSMImageZero);
					} else {
						isc0[dstx] = isr1 + IMUL6(isr0);
					}
					if (ap) {
						ip += isc0[dstx];
						ip /= SKIPSMImagePixelType(ap);
						da.set(DestPixelType(ip), dx);
						daa.set(DestAlphaMax, dax);
					} else {
						da.set(DestImageZero, dx);
						daa.set(DestAlphaZero, dax);
					}
				} else {
					// Previous srcx was odd
					SKIPSMAlphaPixelType ap = asc1[dstx] + AMUL6(asc0[dstx]) + ascp[dstx];
					asc1[dstx] = asc0[dstx] + ascp[dstx];
					if (wraparound) {
						asc0[dstx] = asr1 + AMUL6(asr0) + asrp + (aa(ay) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
					} else {
						asc0[dstx] = asr1 + AMUL6(asr0) + asrp;
					}
					ap += asc0[dstx];

					SKIPSMImagePixelType ip = isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx];
					isc1[dstx] = isc0[dstx] + iscp[dstx];
					if (wraparound) {
						isc0[dstx] = isr1 + IMUL6(isr0) + isrp + (aa(ay) ? SKIPSMImagePixelType(sa(sy)) : SKIPSMImageZero);
					} else {
						isc0[dstx] = isr1 + IMUL6(isr0) + isrp;
					}
					if (ap) {
						ip += isc0[dstx];
						ip /= SKIPSMImagePixelType(ap);
						da.set(DestPixelType(ip), dx);
						daa.set(DestAlphaMax, dax);
					} else {
						da.set(DestImageZero, dx);
						daa.set(DestAlphaZero, dax);
					}
				}

				++dy.y;
				++day.y;
			} else {
				// First entry in odd-numbered row
				sx = sy;
				ax = ay;
				if (wraparound) {
					asr1 = asr0 + asrp;
					isr1 = isr0 + isrp;
				}
				asr0 = aa(ax) ? SKIPSMAlphaOne : SKIPSMAlphaZero;
				isr0 = aa(ax) ? SKIPSMImagePixelType(sa(sx)) : SKIPSMImageZero;
				// isc*[0] are never used
				++sx.x;
				++ax.x;

				// Main entries in odd-numbered row
				for (evenX = false, srcx = 1, dstx = 0; srcx < src_w; ++srcx, ++sx.x, ++ax.x) {
					SKIPSMAlphaPixelType mcurrent(aa(ax) ? SKIPSMAlphaOne : SKIPSMAlphaZero);
					SKIPSMImagePixelType icurrent(aa(ax) ? SKIPSMImagePixelType(sa(sx)) : SKIPSMImageZero);
					if (evenX) {
						ascp[dstx] = (asr1 + AMUL6(asr0) + asrp + mcurrent) * 4;
						asr1 = asr0 + asrp;
						asr0 = mcurrent;
						iscp[dstx] = (isr1 + IMUL6(isr0) + isrp + icurrent) * 4;
						isr1 = isr0 + isrp;
						isr0 = icurrent;
					} else {
						asrp = mcurrent * 4;
						isrp = icurrent * 4;
						++dstx;
					}
					evenX = !evenX;
				}

				// Last entries in row
				if (!evenX) {
					// previous srcx was even
					++dstx;
					if (wraparound) {
						ascp[dstx] = (asr1 + AMUL6(asr0) + (aa(ay) ? (SKIPSMAlphaOne * 4) : SKIPSMAlphaZero)
							+ (aa(ay, Diff2D(1,0)) ? SKIPSMAlphaOne : SKIPSMAlphaZero)
							) * 4;
						iscp[dstx] = (isr1 + IMUL6(isr0) + (aa(ay) ? (SKIPSMImagePixelType(sa(sy)) * 4) : SKIPSMImageZero)
							+ (aa(ay, Diff2D(1, 0)) ? SKIPSMImagePixelType(sa(sy, Diff2D(1, 0))) : SKIPSMImageZero)
							) * 4;
					} else {
						ascp[dstx] = (asr1 + AMUL6(asr0)) * 4;
						iscp[dstx] = (isr1 + IMUL6(isr0)) * 4;
					}
				} else {
					// previous srcx was odd
					if (wraparound) {
						ascp[dstx] = (asr1 + AMUL6(asr0) + asrp + (aa(ay) ? SKIPSMAlphaOne : SKIPSMAlphaZero)) * 4;
						iscp[dstx] = (isr1 + IMUL6(isr0) + isrp + (aa(ay) ? SKIPSMImagePixelType(sa(sy)) : SKIPSMImageZero)) * 4;
					} else {
						ascp[dstx] = (asr1 + AMUL6(asr0) + asrp) * 4;
						iscp[dstx] = (isr1 + IMUL6(isr0) + isrp) * 4;
					}
				}
			}
			evenY = !evenY;
		}
	}

	{
		if (!evenY) {

			for (dstx = 1, dx = dy, dax = day; dstx < dst_w + 1; ++dstx, ++dx.x, ++dax.x) {
				SKIPSMAlphaPixelType ap = asc1[dstx] + AMUL6(asc0[dstx]);
				if (ap) {
					SKIPSMImagePixelType ip = (isc1[dstx] + IMUL6(isc0[dstx])) / SKIPSMImagePixelType(ap);
					da.set(DestPixelType(ip), dx);
					daa.set(DestAlphaMax, dax);
				} else {
					da.set(DestImageZero, dx);
					daa.set(DestAlphaZero, dax);
				}
			}
		} else {

			for (dstx = 1, dx = dy, dax = day; dstx < dst_w + 1; ++dstx, ++dx.x, ++dax.x) {
				SKIPSMAlphaPixelType ap = asc1[dstx] + AMUL6(asc0[dstx]) + ascp[dstx];
				if (ap) {
					SKIPSMImagePixelType ip = (isc1[dstx] + IMUL6(isc0[dstx]) + iscp[dstx]) / SKIPSMImagePixelType(ap);
					da.set(DestPixelType(ip), dx);
					daa.set(DestAlphaMax, dax);
				} else {
					da.set(DestImageZero, dx);
					daa.set(DestAlphaZero, dax);
				}
			}
		}
	}

	delete [] isc0;
	delete [] isc1;
	delete [] iscp;

	delete [] asc0;
	delete [] asc1;
	delete [] ascp;
}


template <typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType,
	typename SrcImageIterator, typename SrcAccessor,
	typename AlphaIterator, typename AlphaAccessor,
	typename DestImageIterator, typename DestAccessor,
	typename DestAlphaIterator, typename DestAlphaAccessor>
	inline void
	reduce(bool wraparound,
	triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
	pair<AlphaIterator, AlphaAccessor> mask,
	triple<DestImageIterator, DestImageIterator, DestAccessor> dest,
	triple<DestAlphaIterator, DestAlphaIterator, DestAlphaAccessor> destMask)
{
	reduce<SKIPSMImagePixelType, SKIPSMAlphaPixelType>(wraparound,
		src.first, src.second, src.third,
		mask.first, mask.second,
		dest.first, dest.second, dest.third,
		destMask.first, destMask.second, destMask.third);
};

template <class SrcImageIterator, class SrcAccessor,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	transformImageMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
	DestImageIterator dest_upperleft, DestAccessor dest_acc,
	const Functor& func)
{
	vigra::transformImage(src_upperleft, src_lowerright, src_acc,
		dest_upperleft, dest_acc,
		func);
}

template <typename SrcPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertScalarToPyramidFunctor
{

public:
	ConvertScalarToPyramidFunctor() { }

	inline PyramidPixelType operator()(const SrcPixelType &v) const
	{
		return doConvert(v, SrcIsIntegral(), PyramidIsIntegral());
	}

protected:

	typedef typename NumericTraits<SrcPixelType>::isIntegral SrcIsIntegral;
	typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

	//将整数类型数据转换成整数的pyramid value类型。
	inline PyramidPixelType doConvert(const SrcPixelType &v, VigraTrueType, VigraTrueType) const 
	{
		return convertIntegerToFixedPoint(v);
	}

	//将整数类型数据转换成实数的pyramid value类型
	inline PyramidPixelType doConvert(const SrcPixelType &v, VigraTrueType, VigraFalseType) const
	{
		return NumericTraits<SrcPixelType>::toRealPromote(v);
	}

	//将实数类型数据转换成整数的pyramid value类型
	inline PyramidPixelType doConvert(const SrcPixelType &v, VigraFalseType, VigraTrueType) const
	{
		return convertDoubleToFixedPoint(v);
	}

	//将实数类型数据转换成实数的pyramid value类型
	inline PyramidPixelType doConvert(const SrcPixelType &v, VigraFalseType, VigraFalseType) const 
	{
		return v;
	}

	inline PyramidPixelType convertDoubleToFixedPoint(const double &v) const 
	{
		return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << PyramidFractionBits));
	};

	inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) const
	{
		return (PyramidPixelType)v << PyramidFractionBits;
	};

};

template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToPyramidFunctor 
{
	typedef typename SrcVectorType::value_type SrcComponentType;
	typedef typename PyramidVectorType::value_type PyramidComponentType;
	typedef ConvertScalarToPyramidFunctor<SrcComponentType, PyramidComponentType,
		PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
	ConvertVectorToPyramidFunctor() : cf() {}

	inline PyramidVectorType operator()(const SrcVectorType &v) const
	{
		return PyramidVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
	}

protected:
	ConvertFunctorType cf;
};
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
void copyToPyramidImage(
	typename SrcImageType::const_traverser src_upperleft,
	typename SrcImageType::const_traverser src_lowerright,
	typename SrcImageType::ConstAccessor sa,
	typename PyramidImageType::traverser dest_upperleft,
	typename PyramidImageType::Accessor da,
	VigraTrueType)
{

	typedef typename SrcImageType::value_type SrcPixelType;
	typedef typename PyramidImageType::value_type PyramidPixelType;

	transformImageMP(src_upperleft, src_lowerright, sa,
		dest_upperleft, da,
		ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType, PyramidIntegerBits, PyramidFractionBits>());
};

template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
void copyToPyramidImage(
	typename SrcImageType::const_traverser src_upperleft,
	typename SrcImageType::const_traverser src_lowerright,
	typename SrcImageType::ConstAccessor sa,
	typename PyramidImageType::traverser dest_upperleft,
	typename PyramidImageType::Accessor da,
	VigraFalseType) {

		typedef typename SrcImageType::value_type SrcVectorType;
		typedef typename PyramidImageType::value_type PyramidVectorType;
/*
		if (UseCIECAM)
		{
			int w = src_lowerright.x - src_upperleft.x;
			int twentyPercent = 1 + ((src_lowerright.y - src_upperleft.y) / 5);
			int tick = 1;
			for (int y = 0; src_upperleft.y < src_lowerright.y; ++src_upperleft.y, ++dest_upperleft.y, ++y) 
			{

				transformLine(src_upperleft.rowIterator(),
					src_upperleft.rowIterator() + w, sa,
					dest_upperleft.rowIterator(), da,
					ConvertVectorToJCHPyramidFunctor<SrcVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
			}

		} */
	
			transformImageMP(src_upperleft, src_lowerright, sa,
				dest_upperleft, da,
				ConvertVectorToPyramidFunctor<SrcVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
		

};

template <typename DestPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToScalarFunctor
{
public:
	ConvertPyramidToScalarFunctor() { }

	inline DestPixelType operator()(const PyramidPixelType &v) const 
	{
		return doConvert(v, DestIsIntegral(), PyramidIsIntegral());
	}

protected:

	typedef typename NumericTraits<DestPixelType>::isIntegral DestIsIntegral;
	typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

	inline DestPixelType doConvert(const PyramidPixelType &v, VigraTrueType, VigraTrueType) const
	{
		PyramidPixelType half = 1 << (PyramidFractionBits-1);
		PyramidPixelType quarter = 1 << (PyramidFractionBits-2);
		PyramidPixelType threeQuarter = 3 << (PyramidFractionBits-2);

		PyramidPixelType vFraction = v & ((1 << PyramidFractionBits) - 1);

		if ((vFraction >= quarter) && (vFraction < threeQuarter)) 
		{
			PyramidPixelType random = (PyramidPixelType(::Twister()) & (half - 1)) + quarter;
			if (random <= vFraction) 
			{
				return DestPixelType(NumericTraits<DestPixelType>::fromPromote((v >> PyramidFractionBits) + 1));
			} else {
				return DestPixelType(NumericTraits<DestPixelType>::fromPromote(v >> PyramidFractionBits));
			}
		}
		else if (vFraction >= quarter)
		{
			return DestPixelType(NumericTraits<DestPixelType>::fromPromote((v >> PyramidFractionBits) + 1));
		}
		else 
		{
			return DestPixelType(NumericTraits<DestPixelType>::fromPromote(v >> PyramidFractionBits));
		}

	}

	//将整数的pyramid value类型数据转换成整数类型。
	inline DestPixelType doConvert(const PyramidPixelType &v, VigraTrueType, VigraFalseType) const 
	{
		double d = dither(v);
		return NumericTraits<DestPixelType>::fromRealPromote(d);
	}
	//将整数的pyramid value类型数据转换成实数类型。
	inline DestPixelType doConvert(const PyramidPixelType &v, VigraFalseType, VigraTrueType) const 
	{
		return convertFixedPointToDouble(v);
	}
	//将实数的pyramid value类型数据转换成实数类型。
	inline DestPixelType doConvert(const PyramidPixelType &v, VigraFalseType, VigraFalseType) const 
	{
		return v;
	}

	inline double dither(const double &v) const 
	{
		double vFraction = v - floor(v);
		if (vFraction > 0.25 && vFraction <= 0.75) 
		{
			double random = 0.5 * (double)::Twister() / UINT_MAX;
			if ((vFraction - 0.25) >= random)
			{
				return ceil(v);
			} 
			else 
			{
				return floor(v);
			}
		}
		else 
		{
			return v;
		}
	}

	inline double convertFixedPointToDouble(const PyramidPixelType &v) const
	{
		return NumericTraits<PyramidPixelType>::toRealPromote(v) / (double)(1 << PyramidFractionBits);
	};

};


template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToVectorFunctor
{

	typedef typename DestVectorType::value_type DestComponentType;
	typedef typename PyramidVectorType::value_type PyramidComponentType;
	typedef ConvertPyramidToScalarFunctor<DestComponentType, PyramidComponentType,
		PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
	ConvertPyramidToVectorFunctor() : cf() {}

	inline DestVectorType operator()(const PyramidVectorType &v) const
	{
		return DestVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
	}

protected:
	ConvertFunctorType cf;
};

template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
	typename PyramidImageType::const_traverser src_upperleft,
	typename PyramidImageType::const_traverser src_lowerright,
	typename PyramidImageType::ConstAccessor sa,
	typename MaskImageType::const_traverser mask_upperleft,
	typename MaskImageType::ConstAccessor ma,
	typename DestImageType::traverser dest_upperleft,
	typename DestImageType::Accessor da,
	VigraTrueType) 
{

	typedef typename DestImageType::value_type DestPixelType;
	typedef typename PyramidImageType::value_type PyramidPixelType;

	transformImageIfMP(src_upperleft, src_lowerright, sa,
		mask_upperleft, ma,
		dest_upperleft, da,
		ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType, PyramidIntegerBits, PyramidFractionBits>());
};


template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
	typename PyramidImageType::const_traverser src_upperleft,
	typename PyramidImageType::const_traverser src_lowerright,
	typename PyramidImageType::ConstAccessor sa,
	typename MaskImageType::const_traverser mask_upperleft,
	typename MaskImageType::ConstAccessor ma,
	typename DestImageType::traverser dest_upperleft,
	typename DestImageType::Accessor da,
	VigraFalseType) {

		typedef typename DestImageType::value_type DestVectorType;
		typedef typename PyramidImageType::value_type PyramidVectorType;



		transformImageIfMP(src_upperleft, src_lowerright, sa,
			mask_upperleft, ma,
			dest_upperleft, da,
			ConvertPyramidToVectorFunctor<DestVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());

};


template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
	typename PyramidImageType::const_traverser src_upperleft,
	typename PyramidImageType::const_traverser src_lowerright,
	typename PyramidImageType::ConstAccessor sa,
	typename MaskImageType::const_traverser mask_upperleft,
	typename MaskImageType::ConstAccessor ma,
	typename DestImageType::traverser dest_upperleft,
	typename DestImageType::Accessor da)
{

	typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

	copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType, PyramidIntegerBits, PyramidFractionBits>(
		src_upperleft,
		src_lowerright,
		sa,
		mask_upperleft,
		ma,
		dest_upperleft,
		da,
		src_is_scalar());
};
template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
	triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
	pair<typename MaskImageType::const_traverser, typename MaskImageType::ConstAccessor> mask,
	pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest) {

		copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType, PyramidIntegerBits, PyramidFractionBits>(
			src.first,
			src.second,
			src.third,
			mask.first,
			mask.second,
			dest.first,
			dest.second);
};



template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyToPyramidImage(
	typename SrcImageType::const_traverser src_upperleft,
	typename SrcImageType::const_traverser src_lowerright,
	typename SrcImageType::ConstAccessor sa,
	typename PyramidImageType::traverser dest_upperleft,
	typename PyramidImageType::Accessor da) {

		typedef typename NumericTraits<typename SrcImageType::value_type>::isScalar src_is_scalar;

		copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>(
			src_upperleft,
			src_lowerright,
			sa,
			dest_upperleft,
			da,
			src_is_scalar());
};
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType,
	int PyramidIntegerBits, int PyramidFractionBits,
	typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType>
	vector<PyramidImageType*>*
	gaussianPyramid(unsigned int numLevels,
	bool wraparound,
	typename SrcImageType::const_traverser src_upperleft,
	typename SrcImageType::const_traverser src_lowerright,
	typename SrcImageType::ConstAccessor sa,
	typename AlphaImageType::const_traverser alpha_upperleft,
	typename AlphaImageType::ConstAccessor aa)
{
	vector<PyramidImageType*>* gp = new vector<PyramidImageType*>();

	int w = src_lowerright.x - src_upperleft.x;
	int h = src_lowerright.y - src_upperleft.y;

	PyramidImageType* gp0 = new PyramidImageType(w, h);

	copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>
		(src_upperleft, src_lowerright, sa, gp0->upperLeft(), gp0->accessor());

	gp->push_back(gp0);


	PyramidImageType* lastGP = gp0;
	AlphaImageType* lastA = NULL;
	for (unsigned int l = 1; l < numLevels; l++)
	{
		w = (w + 1) >> 1;
		h = (h + 1) >> 1;

		PyramidImageType* gpn = new PyramidImageType(w, h);
		AlphaImageType* nextA = new AlphaImageType(w, h);

		if (lastA == NULL)
		{
			reduce<SKIPSMImagePixelType, SKIPSMAlphaPixelType>
				(wraparound,
				srcImageRange(*lastGP), maskIter(alpha_upperleft, aa),
				destImageRange(*gpn), destImageRange(*nextA));
		}
		else 
		{
			reduce<SKIPSMImagePixelType, SKIPSMAlphaPixelType>
				(wraparound,
				srcImageRange(*lastGP), maskImage(*lastA),
				destImageRange(*gpn), destImageRange(*nextA));
		}
		
		gp->push_back(gpn);
		lastGP = gpn;
		delete lastA;
		lastA = nextA;
	}

	delete lastA;

	return gp;
}

template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType,
	int PyramidIntegerBits, int PyramidFractionBits,
	typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType>
	inline vector<PyramidImageType*>*
	gaussianPyramid(unsigned int numLevels,
	bool wraparound,
	triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
	pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha)
{
	return gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits, SKIPSMImagePixelType, SKIPSMAlphaPixelType>
		(numLevels, wraparound,
		src.first, src.second, src.third,
		alpha.first, alpha.second);
}

template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType,
	int PyramidIntegerBits, int PyramidFractionBits,
	typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType>
	vector<PyramidImageType*>*
	laplacianPyramid(const char* exportName, unsigned int numLevels,
	bool wraparound,
	typename SrcImageType::const_traverser src_upperleft,
	typename SrcImageType::const_traverser src_lowerright,
	typename SrcImageType::ConstAccessor sa,
	typename AlphaImageType::const_traverser alpha_upperleft,
	typename AlphaImageType::ConstAccessor aa)
{
	vector <PyramidImageType*>* gp;
	gp =  gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits, SKIPSMImagePixelType, SKIPSMAlphaPixelType>
		(numLevels, wraparound,
		src_upperleft, src_lowerright, sa,
		alpha_upperleft, aa);


	for (unsigned int l = 0; l < (numLevels-1); l++) 
	{
		expand<SKIPSMImagePixelType>(false, wraparound,
			srcImageRange(*((*gp)[l+1])),
			destImageRange(*((*gp)[l])));

	}

	return gp;
}
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType,
	int PyramidIntegerBits, int PyramidFractionBits,
	typename SKIPSMImagePixelType, typename SKIPSMAlphaPixelType>
	inline vector<PyramidImageType*>*
	laplacianPyramid(const char* exportName, unsigned int numLevels,
	bool wraparound,
	triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
	pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha)
{
	return laplacianPyramid<SrcImageType, AlphaImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits, SKIPSMImagePixelType, SKIPSMAlphaPixelType>
		(exportName,
		numLevels, wraparound,
		src.first, src.second, src.third,
		alpha.first, alpha.second);
}


template <typename ImagePixelComponentType>
unsigned int
	filterHalfWidth(const unsigned int levels)
{
	vigra_precondition(levels >= 1 && levels <= MAX_PYRAMID_LEVELS,
		"filterHalfWidth: levels outside of permissible range");

	// This is the arithmetic half width.
	const int halfWidth = levels == 1 ? 0 : (1 << (levels + 1)) - 4;

	return halfWidth;
}

template <typename ImagePixelComponentType>
unsigned int
	roiBounds(const Rect2D &inputUnion,
	const Rect2D &iBB, const Rect2D &mBB, const Rect2D &uBB,
	Rect2D &roiBB,
	bool wraparoundForMask)
{
	unsigned int levels = 1;    
	unsigned int allowableLevels;
	if (ExactLevels <= 0) 
	{
		const unsigned int shortDimension = min(iBB.width(), iBB.height());
		while (levels <= MAX_PYRAMID_LEVELS &&
			(2 * filterHalfWidth<ImagePixelComponentType>(levels) <= shortDimension))
		{
			++levels;
		}

		if (levels == 1)
		{

		}

		if (static_cast<int>(levels) + ExactLevels >= 1) 
		{
			levels += ExactLevels;
		}
		else 
		{
			levels = 1;

		}
	}
	else 
	{
		levels = ExactLevels;
	}

	roiBB = mBB;
	roiBB.addBorder(filterHalfWidth<ImagePixelComponentType>(levels));

	if (wraparoundForMask &&
		(roiBB.left() < 0 || roiBB.right() > uBB.right()))
	{
		roiBB.setUpperLeft(Point2D(0, roiBB.top()));
		roiBB.setLowerRight(Point2D(uBB.right(), roiBB.bottom()));
	}

	roiBB &= uBB;

	unsigned int roiShortDimension = min(roiBB.width(), roiBB.height());
	
	for (allowableLevels = 1; allowableLevels <= levels; allowableLevels++)
	{
		if (roiShortDimension <= 8) 
		{
			break;
		}
		roiShortDimension = (roiShortDimension + 1) >> 1;
	}
	
	return allowableLevels;
}


template <typename SrcIterator,
	typename SrcAccessor,
	typename DesIterator,
	typename DesAccessor>

	void Substract(SrcIterator src1_ul,
	SrcIterator src1_dr,
	SrcAccessor src1_a,
	SrcIterator src2_ul,
	SrcIterator src2_dr,
	SrcAccessor src2_a,
	DesIterator des_ul,
	DesIterator des_dr,
	DesAccessor des_a
	)
{
	int w  = des_dr.x - des_ul.x;
	int h  = des_dr.y - des_ul.y;
#ifdef _OMP_
#pragma omp parallel
	{
#pragma omp for
#endif
		for(int j = 0; j < h; j++)
		{
			for(int i  = 0; i < w; i++)
			{
				SrcIterator src1_curi = src1_ul + Diff2D(i,j);
				SrcIterator src2_curi = src2_ul + Diff2D(i,j);
				//float vv1 = src1_a.(src1_curi);
				//float vv2 = src2_a.(src2_curi);

				DesIterator des_curi = des_ul + Diff2D(i,j);

				des_a.set((src1_a(src1_curi) + src2_a(src2_curi)),des_curi);
			}
		}
#ifdef _OMP_
	}
#endif
}

template <typename SrcIterator,
	typename SrcAccessor,
	typename DesIterator,
	typename DesAccessor>

	void Substract(triple<SrcIterator,SrcIterator,SrcAccessor> src1,
	triple<SrcIterator, SrcIterator, SrcAccessor> src2,
	triple<DesIterator,DesIterator, DesAccessor> dst)
{
	Substract(src1.first,
		src1.second,
		src1.third,
		src2.first,
		src2.second,
		src2.third,
		dst.first,
		dst.second,
		dst.third
		);

}

template <class SrcImageIterator1, class SrcAccessor1,
class SrcImageIterator2, class SrcAccessor2,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	combineTwoImagesMP(SrcImageIterator1 src1_upperleft,
	SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
	SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
	DestImageIterator dest_upperleft, DestAccessor dest_acc,
	const Functor& func)
{
	vigra::combineTwoImages(src1_upperleft, src1_lowerright, src1_acc,
		src2_upperleft, src2_acc,
		dest_upperleft, dest_acc,
		func);
}

template <class SrcImageIterator1, class SrcAccessor1,
class SrcImageIterator2, class SrcAccessor2,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	combineTwoImagesMP(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
	vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
	vigra::pair<DestImageIterator, DestAccessor> dest1,
	const Functor& func)
{
	combineTwoImagesMP(src1.first, src1.second, src1.third,
		src2.first, src2.second,
		dest1.first, dest1.second,
		func);
}


/*
template <class SrcImageIterator1, class SrcAccessor1,
class SrcImageIterator2, class SrcAccessor2,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	combineTwoImagesMP(SrcImageIterator1 src1_upperleft,
	SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
	SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
	DestImageIterator dest_upperleft, DestAccessor dest_acc,
	const Functor& func)
{
	vigra::combineTwoImages(src1_upperleft, src1_lowerright, src1_acc,
		src2_upperleft, src2_acc,
		dest_upperleft, dest_acc,
		func);
}
*/

template <class SrcImageIterator1, class SrcAccessor1,
class SrcImageIterator2, class SrcAccessor2,
class MaskImageIterator, class MaskAccessor,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	combineTwoImagesIfMP(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
	SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
	MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
	DestImageIterator dest_upperleft, DestAccessor dest_acc,
	const Functor& func)
{
	vigra::combineTwoImagesIf(src1_upperleft, src1_lowerright, src1_acc,
		src2_upperleft, src2_acc,
		mask_upperleft, mask_acc,
		dest_upperleft, dest_acc,
		func);
}


template <class SrcImageIterator1, class SrcAccessor1,
class SrcImageIterator2, class SrcAccessor2,
class MaskImageIterator, class MaskAccessor,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	combineTwoImagesIfMP(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
	vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
	vigra::pair<MaskImageIterator, MaskAccessor> mask,
	vigra::pair<DestImageIterator, DestAccessor> dest,
	const Functor& func)
{
	combineTwoImagesIfMP(src1.first, src1.second, src1.third,
		src2.first, src2.second,
		mask.first, mask.second,
		dest.first, dest.second,
		func);
}


static inline double gaussDistribution(double x, double mu, double sigma)
{
	return exp(-0.5 * square((x - mu) / sigma));
}

template <typename InputType, typename ResultType>
class SaturationFunctor 
{
public:
	typedef ResultType result_type;

	SaturationFunctor(double w) : weight(w) {}

	inline ResultType operator()(const InputType& a) const 
	{
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		return f(a, srcIsScalar());
	}

protected:
	template <typename T>
	inline ResultType f(const T& a, VigraTrueType) const 
	{
		return NumericTraits<ResultType>::zero();
	}
	template <typename T>
	inline ResultType f(const T& a, VigraFalseType) const 
	{
		typedef typename T::value_type value_type;
		typedef NumericTraits<value_type> value_traits;
		typedef NumericTraits<ResultType> result_traits;

		const value_type max = std::max(a.red(), std::max(a.green(), a.blue()));
		const value_type min = std::min(a.red(), std::min(a.green(), a.blue()));
		if (max == min)
		{
			return result_traits::zero();
		}
		else
		{
			const double max_value =
				value_traits::isIntegral::asBool ?
				static_cast<double>(value_traits::max()) :
			1.0;
			const double sum = static_cast<double>(max) + static_cast<double>(min);
			const double difference = static_cast<double>(max) - static_cast<double>(min);
			const double saturation =
				sum <= max_value ?
				difference / sum :
			difference / (2.0 * max_value - sum);
			return result_traits::fromRealPromote(weight * saturation);
		}
	}

	double weight;
};

template <class SrcImageIterator, class SrcAccessor,
class MaskImageIterator, class MaskAccessor,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	transformImageIfMP(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
	vigra::pair<MaskImageIterator, MaskAccessor> mask,
	vigra::pair<DestImageIterator, DestAccessor> dest,
	const Functor& func)
{
	transformImageIfMP(src.first, src.second, src.third,
		mask.first, mask.second,
		dest.first, dest.second,
		func);
}

template <class SrcImageIterator, class SrcAccessor,
class MaskImageIterator, class MaskAccessor,
class DestImageIterator, class DestAccessor,
class Functor>
	inline void
	transformImageIfMP(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
	MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
	DestImageIterator dest_upperleft, DestAccessor dest_acc,
	const Functor& func)
{
	vigra::transformImageIf(src_upperleft, src_lowerright, src_acc,
		mask_upperleft, mask_acc,
		dest_upperleft, dest_acc,
		func);
}



template <typename InputType, typename InputAccessor, typename ResultType>
class ExposureFunctor
{
public:
	typedef ResultType result_type;

	ExposureFunctor(double w, double m, double s, InputAccessor a) :
	weight(w), mu(m), sigma(s), acc(a) {}

	inline ResultType operator()(const InputType& a) const 
	{
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		return f(a, srcIsScalar());
	}

protected:

	template <typename T>
	inline ResultType f(const T& a, VigraTrueType) const 
	{
		typedef typename NumericTraits<T>::RealPromote RealType;
		const RealType ra = NumericTraits<T>::toRealPromote(a);
		const double b = NumericTraits<T>::max() * mu;
		const double c = NumericTraits<T>::max() * sigma;
		return NumericTraits<ResultType>::fromRealPromote(weight *gaussDistribution(ra, b, c));
	}


	template <typename T>
	inline ResultType f(const T& a, VigraFalseType) const
	{
		return f(acc.operator()(a), VigraTrueType());
	}

	double weight;
	double mu;
	double sigma;
	InputAccessor acc;
};
template <typename InputType, typename ScaleType, typename ResultType>
class ContrastFunctor 
{
public:
	typedef ResultType result_type;

	ContrastFunctor(double w) : weight(w) {}

	inline ResultType operator()(const InputType& a) const 
	{
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		typedef typename NumericTraits<ScaleType>::isIntegral scaleIsIntegral;
		return f(a, srcIsScalar(), scaleIsIntegral());
	}

protected:
	template <typename T>
	inline ResultType f(const T& a, VigraTrueType, VigraTrueType) const 
	{
		const typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
		return NumericTraits<ResultType>::fromRealPromote(weight * ra / NumericTraits<ScaleType>::max());
	}
	template <typename T>
	inline ResultType f(const T& a, VigraTrueType, VigraFalseType) const 
	{
		const typename NumericTraits<T>::RealPromote ra = NumericTraits<T>::toRealPromote(a);
		return NumericTraits<ResultType>::fromRealPromote(weight * ra);
	}

	template <typename T>
	inline ResultType f(const T& a, VigraFalseType, VigraTrueType) const
	{
		typedef typename T::value_type TComponentType;
		typedef typename NumericTraits<TComponentType>::RealPromote RealTComponentType;
		typedef typename ScaleType::value_type ScaleComponentType;
		const RealTComponentType ra = static_cast<RealTComponentType>(a.lightness());
		return NumericTraits<ResultType>::fromRealPromote(weight * ra / NumericTraits<ScaleComponentType>::max());
	}

	template <typename T>
	inline ResultType f(const T& a, VigraFalseType, VigraFalseType) const 
	{
		typedef typename T::value_type TComponentType;
		typedef typename NumericTraits<TComponentType>::RealPromote RealTComponentType;
		const RealTComponentType ra = static_cast<RealTComponentType>(a.lightness());
		return NumericTraits<ResultType>::fromRealPromote(weight * ra);
	}

	double weight;
};

template <typename InputType, typename ResultType>
class MultiGrayscaleAccessor
{
public:
	typedef ResultType value_type;

	MultiGrayscaleAccessor(const std::string& accessorName) {
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		initializeTypeSpecific(srcIsScalar());
		initialize(accessorName);
	}

	ResultType operator()(const InputType& x) const {
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		return f(x, srcIsScalar());
	}

	template <class Iterator>
	ResultType operator()(const Iterator& i) const {
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		return f(i, srcIsScalar());
	}

	template <class Iterator, class Difference>
	ResultType operator()(const Iterator& i, Difference d) const {
		typedef typename NumericTraits<InputType>::isScalar srcIsScalar;
		return f(i, d, srcIsScalar());
	}

	static const std::string defaultGrayscaleAccessorName() {return "average";}

private:
	typedef enum {AVERAGE, LSTAR, LIGHTNESS, VALUE, LUMINANCE, MIXER} AccKindType;
	typedef std::map<std::string, AccKindType> NameMapType;
	typedef typename NameMapType::const_iterator NameMapConstIterType;

#define CHANNEL_MIXER "channel-mixer"

	void initializeAccessorNameMap() {
		nameMap["average"] = AVERAGE;
		nameMap["l-star"] = LSTAR;
		nameMap["lightness"] = LIGHTNESS;
		nameMap["value"] = VALUE;
		nameMap["luminance"] = LUMINANCE;
		nameMap[CHANNEL_MIXER] = MIXER;
	}

	void initialize(const std::string& accessorName) {
		initializeAccessorNameMap();
		if (accessorName.empty())
		{
			kind = nameMap[defaultGrayscaleAccessorName()];
		}
		else
		{
			std::string name(accessorName);
#ifdef   WIN32
			std::transform(name.begin(), name.end(), name.begin(), tolower);
#else
			//stlp_std::transform(name.begin(), name.end(), name.begin(), tolower);
#endif
			NameMapConstIterType const k = nameMap.find(name);
			if (k == nameMap.end())
			{
				char dummy;
				double red, green, blue;
				if (sscanf(name.c_str(),
					CHANNEL_MIXER "%["
					NUMERIC_OPTION_DELIMITERS "]%lf%["
					NUMERIC_OPTION_DELIMITERS "]%lf%["
					NUMERIC_OPTION_DELIMITERS "]%lf",
					&dummy, &red, &dummy, &green, &dummy, &blue) == 6)
				{
					check_weights(red, green, blue);
					const double sum = red + green + blue;
					redWeight = red / sum;
					greenWeight = green / sum;
					blueWeight = blue / sum;
					kind = MIXER;
				}
				else
				{

					exit(1);
				}

			}
			else
			{
				kind = k->second;
				if (kind == MIXER)
				{

					exit(1);
				}
			}
		}
	}

	void check_weights(double red, double green, double blue) const {
		// TODO: check for isnormal(WEIGHT) before comparison
		if (red < 0.0)
		{

			exit(1);
		}
		if (green < 0.0)
		{

			exit(1);
		}
		if (blue < 0.0)
		{

			exit(1);
		}
		if (red + green + blue == 0.0)
		{

			exit(1);
		}
	}

	void initializeTypeSpecific(VigraTrueType) {}

	void initializeTypeSpecific(VigraFalseType) {
		typedef typename InputType::value_type ValueType;
		labfun = vigra::RGB2LabFunctor<double>(NumericTraits<ValueType>::max());
	}

	ResultType project(const InputType& x) const {
		typedef typename InputType::value_type ValueType;
		switch (kind)
		{
		case AVERAGE:
			return NumericTraits<ResultType>::fromRealPromote
				((NumericTraits<ValueType>::toRealPromote(x.red()) +
				NumericTraits<ValueType>::toRealPromote(x.green()) +
				NumericTraits<ValueType>::toRealPromote(x.blue())) /
				3.0);
		case LSTAR:
			{
				typedef typename vigra::RGB2LabFunctor<double>::result_type LABResultType;
				const LABResultType y = labfun.operator()(x) / 100.0;
				return NumericTraits<ResultType>::fromRealPromote(NumericTraits<ValueType>::max() * y[0]);
			}
		case LIGHTNESS:
			return NumericTraits<ResultType>::fromRealPromote
				((std::min(x.red(), std::min(x.green(), x.blue())) +
				std::max(x.red(), std::max(x.green(), x.blue()))) /
				2.0);
		case VALUE:
			return std::max(x.red(), std::max(x.green(), x.blue()));
		case LUMINANCE:
			return NumericTraits<ResultType>::fromRealPromote(x.luminance());
		case MIXER:
			return NumericTraits<ResultType>::fromRealPromote
				(redWeight * NumericTraits<ValueType>::toRealPromote(x.red()) +
				greenWeight * NumericTraits<ValueType>::toRealPromote(x.green()) +
				blueWeight * NumericTraits<ValueType>::toRealPromote(x.blue()));
		}

		// never reached
		return ResultType();
	}

	// RGB
	ResultType f(const InputType& x, VigraFalseType) const {
		return project(x);
	}

	template <class Iterator>
	ResultType f(const Iterator& i, VigraFalseType) const {
		return project(*i);
	}

	template <class Iterator, class Difference>
	ResultType f(const Iterator& i, Difference d, VigraFalseType) const {
		return project(i[d]);
	}

	// grayscale
	ResultType f(const InputType& x, VigraTrueType) const {return x;}

	template <class Iterator>
	ResultType f(const Iterator& i, VigraTrueType) const {return *i;}

	template <class Iterator, class Difference>
	ResultType f(const Iterator& i, Difference d, VigraTrueType) const {return i[d];}

	NameMapType nameMap;
	AccKindType kind;
	double redWeight, greenWeight, blueWeight;
	vigra::RGB2LabFunctor<double> labfun;
};

/********************************************************/
/*                                                      */
/*                       initLine                       */
/*                                                      */
/********************************************************/

template <class DestIterator, class DestAccessor, class VALUETYPE>
void
	initLineImpl(DestIterator d, DestIterator dend, DestAccessor dest,
	VALUETYPE v, VigraFalseType)
{
	for(; d != dend; ++d)
		dest.set(v, d);
}

template <class DestIterator, class DestAccessor, class FUNCTOR>
void
	initLineImpl(DestIterator d, DestIterator dend, DestAccessor dest,
	FUNCTOR const & f, VigraTrueType)
{
	for(; d != dend; ++d)
		dest.set(f(), d);
}

template <class DestIterator, class DestAccessor, class VALUETYPE>
inline void
	initLine(DestIterator d, DestIterator dend, DestAccessor dest,
	VALUETYPE v)
{
	initLineImpl(d, dend, dest, v, typename FunctorTraits<VALUETYPE>::isInitializer());
}


template <class ImageIterator, class Accessor, class VALUETYPE>
void
	initImage(ImageIterator upperleft, ImageIterator lowerright, 
	Accessor a,  VALUETYPE v)
{
	int w = lowerright.x - upperleft.x;

	for(; upperleft.y < lowerright.y; ++upperleft.y)
	{
		initLine(upperleft.rowIterator(), 
			upperleft.rowIterator() + w, a, v);
	}
}

template <class ImageIterator, class Accessor, class VALUETYPE>
inline 
	void
	initImage(triple<ImageIterator, ImageIterator, Accessor> img, VALUETYPE v)
{
	initImage(img.first, img.second, img.third, v);
}

template <typename DestIterator, typename DestAccessor,
	typename AlphaIterator, typename AlphaAccessor>
	void
	import(const ImageImportInfo& info,
	const pair<DestIterator, DestAccessor>& image,
	const pair<AlphaIterator, AlphaAccessor>& alpha)
{
	typedef typename DestIterator::value_type ImagePixelType;
	typedef typename hdrNumericTraits<ImagePixelType>::PixelType
		ImagePixelComponentType;
	typedef typename AlphaIterator::value_type AlphaPixelType;

	const Diff2D extent = Diff2D(info.width(), info.height());
	const std::string pixelType = info.getPixelType();
	const range_t inputRange = rangeOfPixelType(pixelType);

	/*if (info.numExtraBands() > 0)
	{     
		WriteFunctorAccessor<Threshold<ImagePixelComponentType, AlphaPixelType>, AlphaAccessor>
			ata(Threshold<ImagePixelComponentType, AlphaPixelType>
			(inputRange.second / 2,
			inputRange.second,
			AlphaTraits<AlphaPixelType>::zero(),
			AlphaTraits<AlphaPixelType>::max()),
			alpha.second);
			
		//importImageAlpha(info, image, destIter(alpha.first, ata));
	} */
	importImage(info, image.first, image.second);
	initImage(srcIterRange(alpha.first, alpha.first + extent, alpha.second),
		AlphaTraits<AlphaPixelType>::max());

	//SrcImageIterator::value_type  max_value = NumericTraits<SrcImageIterator::value_type>::max();
	//SrcImageIterator::value_type  min_value = NumericTraits<SrcImageIterator::value_type>::min();
	const double min = 0.0;
	const double max = 1.0;


}

template <typename ImageType, typename AlphaType>
pair<ImageType*, AlphaType*>
	assemble(list<ImageImportInfo*>& imageInfoList, Rect2D& inputUnion, Rect2D& bb)
{
	typedef typename AlphaType::traverser AlphaIteratorType;
	typedef typename AlphaType::Accessor AlphaAccessor;
	if (imageInfoList.empty()) 
	{
		return pair<ImageType*, AlphaType*>(static_cast<ImageType*>(NULL),
			static_cast<AlphaType*>(NULL));
	}
	ImageType* image = new ImageType(inputUnion.size());
	AlphaType* imageA = new AlphaType(inputUnion.size());
	const Diff2D imagePos = imageInfoList.front()->getPosition();

	import(*imageInfoList.front(),
		destIter(image->upperLeft() + imagePos - inputUnion.upperLeft()),
		destIter(imageA->upperLeft() + imagePos - inputUnion.upperLeft()));
	imageInfoList.erase(imageInfoList.begin());
	if (!OneAtATime) 
	{
		list<list<ImageImportInfo*>::iterator> toBeRemoved;

		list<ImageImportInfo*>::iterator i;
		for (i = imageInfoList.begin(); i != imageInfoList.end(); i++) 
		{
			ImageImportInfo* info = *i;

			ImageType* src = new ImageType(info->size());
			AlphaType* srcA = new AlphaType(info->size());

			import(*info, destImage(*src), destImage(*srcA));

			bool overlapFound = false;
			AlphaIteratorType dy =
				imageA->upperLeft() - inputUnion.upperLeft() + info->getPosition();
			AlphaAccessor da = imageA->accessor();
			AlphaIteratorType sy = srcA->upperLeft();
			AlphaIteratorType send = srcA->lowerRight();
			AlphaAccessor sa = srcA->accessor();
			/*for(; sy.y < send.y; ++sy.y, ++dy.y)
			{
				AlphaIteratorType sx = sy;
				AlphaIteratorType dx = dy;
				for(; sx.x < send.x; ++sx.x, ++dx.x) 
				{
					if (sa(sx) && da(dx)) 
					{
						overlapFound = true;
						break;
					}
				}
				if (overlapFound) break;
			}*/

			/*if (!overlapFound) 
			{

				const Diff2D srcPos = info->getPosition();
				copyImageIf(srcImageRange(*src),
					maskImage(*srcA),
					destIter(image->upperLeft() - inputUnion.upperLeft() + srcPos));
				copyImageIf(srcImageRange(*srcA),
					maskImage(*srcA),
					destIter(imageA->upperLeft() - inputUnion.upperLeft() + srcPos));
				toBeRemoved.push_back(i);
			}*/

			delete src;
			delete srcA;
		}

		/*for (list<list<ImageImportInfo*>::iterator>::iterator r = toBeRemoved.begin();
			r != toBeRemoved.end();
			++r) 
		{
			imageInfoList.erase(*r);
		}*/
	}

	FindBoundingRectangle unionRect;
		inspectImageIf(srcIterRange(Diff2D(), Diff2D() + image->size()),
			srcImage(*imageA), unionRect);
		bb = unionRect();
		
	
	return pair<ImageType*, AlphaType*>(image, imageA);
}

template <typename ImageType, typename AlphaType, typename MaskType>
void enfuseMask(triple<typename ImageType::const_traverser, typename ImageType::const_traverser, typename ImageType::ConstAccessor> src,
	pair<typename AlphaType::const_traverser, typename AlphaType::ConstAccessor> mask,
	pair<typename MaskType::traverser, typename MaskType::Accessor> result) {
		typedef typename ImageType::value_type ImageValueType;
		typedef typename ImageType::PixelType PixelType;
		typedef typename NumericTraits<PixelType>::ValueType ScalarType;
		typedef typename MaskType::value_type MaskValueType;

		const typename ImageType::difference_type imageSize = src.second - src.first;
		
		if (WExposure > 0.0)
		{
			typedef MultiGrayscaleAccessor<ImageValueType, ScalarType> MultiGrayAcc;
			MultiGrayAcc ga(GrayscaleProjector);
			ExposureFunctor<ImageValueType, MultiGrayAcc, MaskValueType> ef(WExposure, WMu, WSigma, ga);
			transformImageIfMP(src, mask, result, ef);
		}
	/*	
		if (WContrast > 0.0)
		{
			typedef typename NumericTraits<ScalarType>::Promote LongScalarType;
			typedef IMAGE_TYPE<LongScalarType> GradImage;
			typedef typename GradImage::iterator GradIterator;

			GradImage grad(imageSize);
			MultiGrayscaleAccessor<PixelType, LongScalarType> ga(GrayscaleProjector);
/*
			if (FilterConfig.edgeScale > 0.0)
			{
				GradImage laplacian(imageSize);

				if (FilterConfig.lceScale > 0.0)
				{

					GradImage lce(imageSize);
					gaussianSharpening(src.first, src.second, ga,
						lce.upperLeft(), lce.accessor(),
						FilterConfig.lceFactor, FilterConfig.lceScale);
					laplacianOfGaussian(lce.upperLeft(), lce.lowerRight(), lce.accessor(),
						laplacian.upperLeft(), MagnitudeAccessor<LongScalarType>(),
						FilterConfig.edgeScale);
				}
				else
				{
					laplacianOfGaussian(src.first, src.second, ga,
						laplacian.upperLeft(), MagnitudeAccessor<LongScalarType>(),
						FilterConfig.edgeScale);
				}


				const double minCurve =
					MinCurvature.isPercentage ?
					static_cast<double>(NumericTraits<ScalarType>::max()) * MinCurvature.value / 100.0 :
				MinCurvature.value;
				if (minCurve <= 0.0)
				{

					transformImageIfMP(laplacian.upperLeft(), laplacian.lowerRight(), laplacian.accessor(),
						mask.first, mask.second,
						grad.upperLeft(), grad.accessor(),
						ClampingFunctor<LongScalarType, LongScalarType>
						(static_cast<LongScalarType>(-minCurve), LongScalarType(),
						NumericTraits<LongScalarType>::max(), NumericTraits<LongScalarType>::max()));
				}
				else
				{

					GradImage localContrast(imageSize);

					localStdDevIf(src.first, src.second, ga,
						mask.first, mask.second,
						localContrast.upperLeft(), localContrast.accessor(),
						Size2D(ContrastWindowSize, ContrastWindowSize));

					combineTwoImagesIfMP(laplacian.upperLeft(), laplacian.lowerRight(), laplacian.accessor(),
						localContrast.upperLeft(), localContrast.accessor(),
						mask.first, mask.second,
						grad.upperLeft(), grad.accessor(),
						FillInFunctor<LongScalarType, LongScalarType>
						(static_cast<LongScalarType>(minCurve), 
						1.0, 
						minCurve / NumericTraits<ScalarType>::max())); 
				}
			}
			else
			{

				localStdDevIf(src.first, src.second, ga,
					mask.first, mask.second,
					grad.upperLeft(), grad.accessor(),
					Size2D(ContrastWindowSize, ContrastWindowSize));
			}
			*/
	/*		ContrastFunctor<LongScalarType, ScalarType, MaskValueType> cf(WContrast);
			//combineTwoImagesIfMP(srcImageRange(grad), result, mask, result,
				//const_parameters(bind(cf, boost::lambda::_1) + boost::lambda::_2)); 
			printf("nothing");
		}
*/
		if (WSaturation > 0.0)
		{
			SaturationFunctor<ImageValueType, MaskValueType> sf(WSaturation);
			combineTwoImagesIfMP(src, result, mask, result,
				const_parameters(boost::lambda::bind(sf, boost::lambda::_1) + boost::lambda::_2));

		}
/*
		if (WEntropy > 0.0) 
		{
			typedef typename ImageType::PixelType PixelType;
			typedef typename NumericTraits<PixelType>::ValueType ScalarType;
			typedef IMAGE_TYPE<PixelType> Image;
			Image entropy(imageSize);

			if (EntropyLowerCutoff.value > 0.0 ||
				(EntropyLowerCutoff.isPercentage && EntropyLowerCutoff.value < 100.0) ||
				(!EntropyLowerCutoff.isPercentage && EntropyUpperCutoff.value < NumericTraits<ScalarType>::max()))
			{
				const ScalarType lowerCutoff =
					EntropyLowerCutoff.isPercentage ?
					EntropyLowerCutoff.value *
					static_cast<double>(NumericTraits<ScalarType>::max()) /
					100.0 :
				EntropyLowerCutoff.value;
				const ScalarType upperCutoff =
					EntropyUpperCutoff.isPercentage ?
					EntropyUpperCutoff.value *
					static_cast<double>(NumericTraits<ScalarType>::max()) /
					100.0 :
				EntropyUpperCutoff.value;

				Image trunc(imageSize);
				ClampingFunctor<PixelType, PixelType>
					cf((PixelType(lowerCutoff)),  
					(PixelType(ScalarType())), 
					(PixelType(upperCutoff)),
					(PixelType(NumericTraits<ScalarType>::max())));
				transformImageMP(src.first, src.second, src.third,
					trunc.upperLeft(), trunc.accessor(),
					cf);
				localEntropyIf(trunc.upperLeft(), trunc.lowerRight(), trunc.accessor(),
					mask.first, mask.second,
					entropy.upperLeft(), entropy.accessor(),
					Size2D(EntropyWindowSize, EntropyWindowSize));
			}
			else
			{
				localEntropyIf(src.first, src.second, src.third,
					mask.first, mask.second,
					entropy.upperLeft(), entropy.accessor(),
					Size2D(EntropyWindowSize, EntropyWindowSize));
			}

			EntropyFunctor<PixelType, MaskValueType> ef(WEntropy);
			combineTwoImagesIfMP(srcImageRange(entropy), result, mask, result,
				const_parameters(bind(ef, boost::lambda::_1) + boost::lambda::_2));
		}*/
};


template <typename ImagePixelType>
int hdrMain(const FileNameList& anInputFileNameList,const list<ImageImportInfo*>& anImageInfoList, Rect2D& anInputUnion)
{
	const int imgListInfoProLength = 42;
	const int imgListProLength = 45;
	const int imgLeftProLength = 8;
	printf("Main program is started.\n");
	typedef typename hdrNumericTraits<ImagePixelType>::PixelType ImagePixelComponentType;
	typedef IMAGE_TYPE<ImagePixelType> ImageType;
	typedef typename hdrNumericTraits<ImagePixelType>::ImagePyramidPixelType ImagePyramidPixelType;
	typedef typename hdrNumericTraits<ImagePixelType>::ImagePyramidType ImagePyramidType;
	typedef BasicImage<vigra::UInt8> AlphaType; 
	typedef BasicImage<float> MaskType;
	typedef UInt8 MaskPixelType;
	enum {ImagePyramidIntegerBits = hdrNumericTraits<ImagePixelType>::ImagePyramidIntegerBits};
	enum {ImagePyramidFractionBits = hdrNumericTraits<ImagePixelType>::ImagePyramidFractionBits};
	enum {MaskPyramidIntegerBits = hdrNumericTraits<ImagePixelType>::MaskPyramidIntegerBits};
	enum {MaskPyramidFractionBits = hdrNumericTraits<ImagePixelType>::MaskPyramidFractionBits};
	typedef vigra::Int32 SKIPSMIMAGE ;
	typedef vigra::Int32 SKIPSMMaskPixelType ;
	typedef RGBValue<SKIPSMIMAGE,0,1,2> SKIPSMImagePixelType;
	typedef vigra::Int16 SKIPSMAlphaPixelType;
	typedef typename hdrNumericTraits<ImagePixelType>::MaskPyramidPixelType MaskPyramidPixelType;
	typedef typename hdrNumericTraits<ImagePixelType>::MaskPyramidType MaskPyramidType;
	
#ifdef WIN32

	DWORD  dwGTCBegin, dwGTCEnd;  

#endif

#ifdef WIN32

	dwGTCBegin = GetTickCount();  

#endif
	
	list <triple<ImageType*, AlphaType*, MaskType*> > imageList;

	MaskType *normImage = new MaskType(anInputUnion.size());

	// Result image. Alpha will be union of all input alphas.
	pair<ImageType*, AlphaType*> outputPair(static_cast<ImageType*>(NULL),
		new AlphaType(anInputUnion.size()));
	list<ImageImportInfo*> imageInfoList(anImageInfoList);
	const unsigned numberOfImages = imageInfoList.size();

	int count=0;
	float percent=0.0;
	int w=imageInfoList.size ();


	unsigned m = 0;
	FileNameList::const_iterator inputFileNameIterator(anInputFileNameList.begin());

	int xf_i = 1;
	while (!imageInfoList.empty()) 
	{
		printf("Calculating the coefficient mask for image%d\n",(m+1));
		Rect2D imageBB;
		pair<ImageType*, AlphaType*> imagePair =
			assemble<ImageType, AlphaType>(imageInfoList, anInputUnion, imageBB);
		
		MaskType* mask = new MaskType(anInputUnion.size());
		
		enfuseMask<ImageType, AlphaType, MaskType>(srcImageRange(*(imagePair.first)),
			srcImage(*(imagePair.second)),
			destImage(*mask));
		//exportImage(srcImageRange(*mask),vigra::ImageExportInfo("mask_image.bmp"));

		copyImageIf(srcImageRange(*(imagePair.second)),
			maskImage(*(imagePair.second)),
			destImage(*(outputPair.second)));
		/*combineTwoImagesMP(srcImageRange(*mask),
			srcImage(*normImage),
			destImage(*normImage),
			Arg1() + Arg2());*/
		Substract(srcImageRange(*mask),srcImageRange(*normImage), destImageRange(*normImage));
		/*combineTwoImagesMP(srcImageRange(*mask),
			srcImage(*normImage),
			destImage(*normImage),
			vigra::functor::UnaryFunctor< vigra::functor_add<vigra::functor::UnaryFunctor<vigra::functor::ArgumentFunctor1>,vigra::functor::UnaryFunctor<vigra::functor::ArgumentFunctor2>>>);
		*/
		//exportImage(srcImageRange(*normImage),vigra::ImageExportInfo("norm_image.tif"));
		imageList.push_back(make_triple(imagePair.first, imagePair.second, mask));

		++m;
		++inputFileNameIterator;
/*

		if(pl)
		{
			if(!pl->Update(begin+xf_i*imgListInfoProLength/w))
			{
				return 0;
			}
		}
		++xf_i;

*/
	}
	const int totalImages = imageList.size();
	//typename MaskPixelType maxMaskPixelType ;
	//typename ImagePixelType::value_type  maxMaskPixelType = NumericTraits<ImagePixelType::value_type>::max();
	typename AlphaType::value_type  maxMaskPixelType = NumericTraits<AlphaType::value_type>::max();
		//NumericTraits<typename NumericTraits<ImagePixelType>::MaskPixelType>::max();
	if (UseHardMask)
	{
		const Size2D sz = normImage->size();
		typename list< triple<ImageType*, AlphaType*, MaskType*> >::iterator imageIter;
#ifdef OPENMP
#pragma omp parallel for private (imageIter)
#endif
		for (int y = 0; y < sz.y; ++y)
		{
			for (int x = 0; x < sz.x; ++x) 
			{
				float max = 0.0f;
				int maxi = 0;
				int i = 0;
				for (imageIter = imageList.begin();
					imageIter != imageList.end();
					++imageIter)
				{
					const float w = static_cast<float>((*imageIter->third)(x, y));
					if (w > max)
					{
						max = w;
						maxi = i;
					}
					i++;
				}
				i = 0;
				for (imageIter = imageList.begin();
					imageIter != imageList.end();
					++imageIter) 
				{

					if (max == 0.0f)
					{
						(*imageIter->third)(x, y) =
							static_cast<MaskPixelType>(maxMaskPixelType) / totalImages;
					}
					else if (i == maxi) 
					{
						(*imageIter->third)(x, y) = maxMaskPixelType;
					} 
					else 
					{
						(*imageIter->third)(x, y) = 0.0f;
					}
					
					i++;
				}
			}
		}
		unsigned i = 0;
	}
	Rect2D junkBB;
	const unsigned int numLevels =
		roiBounds<ImagePixelComponentType>(anInputUnion, anInputUnion, anInputUnion, anInputUnion,
		junkBB,
		WrapAround != OpenBoundaries);

	vector<ImagePyramidType*> *resultLP = NULL;

	m = 0;
	int xf_j = 1;
	while (!imageList.empty())
	{
		triple<ImageType*, AlphaType*, MaskType*> imageTriple = imageList.front();
		imageList.erase(imageList.begin());

		//begin = (int)(0.5+percent);
		std::ostringstream oss0;
		oss0 << "imageGP" << m << "_";
		printf("Consruct the Laplasian pyramid for image%d\n",(m+1));
		vector<ImagePyramidType*> *imageLP =
			laplacianPyramid<ImageType, AlphaType, ImagePyramidType,
			ImagePyramidIntegerBits, ImagePyramidFractionBits,
			SKIPSMImagePixelType, SKIPSMAlphaPixelType>(
			oss0.str().c_str(),
			numLevels, WrapAround != OpenBoundaries,
			srcImageRange(*(imageTriple.first)),
			maskImage(*(imageTriple.second)));
		
delete imageTriple.first;
		delete imageTriple.second;
		printf("Consruct the Gaussian pyramid for image%d\n",(m+1));
		vector<MaskPyramidType*> *maskGP;
		maskGP =
			gaussianPyramid<MaskType, AlphaType, MaskPyramidType,
			MaskPyramidIntegerBits, MaskPyramidFractionBits,
			SKIPSMMaskPixelType, SKIPSMAlphaPixelType>
			(numLevels,
			WrapAround != OpenBoundaries,
			srcImageRange(*(imageTriple.third)),
			maskImage(*(outputPair.second)));
			
		delete imageTriple.third;


		ConvertScalarToPyramidFunctor<MaskPixelType,
			MaskPyramidPixelType,
			MaskPyramidIntegerBits,
			MaskPyramidFractionBits> maskConvertFunctor;
		MaskPyramidPixelType maxMaskPyramidPixelValue = maskConvertFunctor(maxMaskPixelType);
		
		for (unsigned int i = 0; i < maskGP->size(); ++i) 
		{
			//exportImage(srcImageRange(*((*imageLP)[i])), vigra::ImageExportInfo("laplacian.bmp"));
			//exportImage(srcImageRange(*((*maskGP)[i])), vigra::ImageExportInfo("gaussian.bmp"));
			printf("combining the the %dth level of laplasian and gaussian pyramid\n",(i+1));
			combineTwoImagesMP(srcImageRange(*((*imageLP)[i])),
				srcImage(*((*maskGP)[i])),
				destImage(*((*imageLP)[i])),
				ImageMaskMultiplyFunctor<MaskPyramidPixelType>(maxMaskPyramidPixelValue));

			
			delete (*maskGP)[i];
		}
		delete maskGP;

		printf("Adding imageLP to resultLP to obtain one result pyramid for image%d\n",(m+1));
		
		if (resultLP != NULL)
		{
			// Add imageLP to resultLP.
			for (unsigned int i = 0; i < imageLP->size(); ++i)
			{
				
				//int w = (*imageLP)[i].width();
				//int h = (*imageLP)[i].height();
				/*combineTwoImagesMP(srcImageRange(*((*imageLP)[i])),
					srcImage(*((*resultLP)[i])),
					destImage(*((*resultLP)[i])),
					Arg1() + Arg2());*/
				//exportImage(srcImageRange(*((*imageLP)[i])), vigra::ImageExportInfo("laplacian.bmp"));
				
				Substract(srcImageRange(*((*imageLP)[i])),srcImageRange(*((*resultLP)[i])), destImageRange(*((*resultLP)[i])));
				//exportImage(srcImageRange(*((*resultLP)[i])), vigra::ImageExportInfo("result.bmp"));
				delete (*imageLP)[i];
				
			}
			delete imageLP;
		}
		else 
		{
			resultLP = imageLP;
		}

		++m;


	
	}

	delete normImage;
	printf("Destructing the combined pyramid and obtaining result enfused image \n");
	collapsePyramid<SKIPSMImagePixelType>(WrapAround != OpenBoundaries, resultLP);

	outputPair.first = new ImageType(anInputUnion.size());

	copyFromPyramidImageIf<ImagePyramidType, AlphaType, ImageType,
		ImagePyramidIntegerBits, ImagePyramidFractionBits>(
		srcImageRange(*((*resultLP)[0])),
		maskImage(*(outputPair.second)),
		destImage(*(outputPair.first)));
	
	for (unsigned int i = 0; i < resultLP->size(); ++i) 
	{
		delete (*resultLP)[i];
	}
	delete resultLP;
	
#ifdef WIN32

	dwGTCEnd = GetTickCount();  
	printf("The image enfusion is don successfully\n");
	//printf("\t\t [timing] ---  cost %d ms\n", dwGTCEnd - dwGTCBegin);  
	FILE* search_file = fopen( "process_time.txt", "a" ); 
	fprintf(search_file,"HDR process time is %d\n",(dwGTCEnd - dwGTCBegin));
	fclose(search_file);//turgun

#endif
	
	exportImage(srcImageRange(*(outputPair.first)),vigra::ImageExportInfo("result.tif"));

	return 0;
}