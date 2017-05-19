#ifndef DEFIER_GRAPH_H
#define DEFIER_GRAPH_H

#include <assert.h>

#include <io.h>
#include <vigra/diff2d.hxx>
#include <vigra/rgbvalue.hxx>


#ifndef _WIN32
#include <unistd.h>
#include <sys/mman.h>
#else
#include <windows.h>
#undef max
#undef min
#endif

namespace vigra{

template <class PIXELTYPE>
class FileMapImageSequentialAccessIteratorPolicy;

template <class IMAGEITERATOR, class IMAGETYPE, class PIXELTYPE, class REFERENCE, class POINTER>
class FileMapImageIteratorBase
{
public:
	typedef FileMapImageIteratorBase<IMAGEITERATOR, IMAGETYPE, PIXELTYPE, REFERENCE, POINTER> self_type;
	typedef IMAGETYPE image_type;
	typedef PIXELTYPE value_type;
	typedef PIXELTYPE PixelType;
	typedef REFERENCE reference;
	typedef REFERENCE index_reference;
	typedef POINTER pointer;
	typedef Diff2D difference_type;
	typedef image_traverser_tag iterator_category;
	typedef RowIterator<IMAGEITERATOR> row_iterator;
    typedef ColumnIterator<IMAGEITERATOR> column_iterator;

	typedef int MoveX;
	typedef int MoveY;
	friend class FileMapImageSequentialAccessIteratorPolicy<IMAGEITERATOR>;
	MoveX x;
	MoveY y;

	IMAGEITERATOR & operator+=(const difference_type & s)
	{
		x += s.x;
		y += s.y;
		return static_cast<IMAGEITERATOR &>(*this);
	}

	IMAGEITERATOR & operator-=(const difference_type & s)
	{
		x -= s.x;
		y -= s.y;
		return static_cast<IMAGEITERATOR &>(*this);
	}

	IMAGEITERATOR operator+(const difference_type & s) const
	{
        IMAGEITERATOR ret(static_cast<const IMAGEITERATOR &>(*this));
        ret += s;
        return ret;
    }

    IMAGEITERATOR operator-(const difference_type & s) const
	{
        IMAGEITERATOR ret(static_cast<const IMAGEITERATOR &>(*this));
        ret -= s;
        return ret;
    }

	difference_type operator-(const FileMapImageIteratorBase & rhs) const
	{
		return difference_type(x-rhs.x, y-rhs.y);
	}

	bool operator==(const FileMapImageIteratorBase & rhs) const
	{
		return (x == rhs.x) && (y == rhs.y);
	}

	bool operator!=(const FileMapImageIteratorBase & rhs) const
	{
		return (x != rhs.x) || (y != rhs.y);
	}

	reference operator*() const
	{
		return (*m_pImage)(x,y);
	}

	pointer operator->() const {
		return (*m_pImage)(x,y);
	}

	index_reference operator[](const difference_type & d) const
	{
		return (*m_pImage)(x+d.x,y+d.y);
	}

	index_reference operator()(int dx, int dy) const
	{
		return (*m_pImage)(x+dx,y+dy);
	}

	pointer operator[](int dy) const
	{
		return (*m_pImage)(x,y+dy);
	}

	row_iterator rowIterator() const
	{
		return row_iterator(static_cast<const IMAGEITERATOR &>(*this));
	}

	column_iterator columnIterator() const
	{
		return column_iterator(static_cast<const IMAGEITERATOR &>(*this));
	}
protected:
	FileMapImageIteratorBase(int X, int Y, image_type * I) : x(X), y(Y), m_pImage(I)
	{
    }

    FileMapImageIteratorBase(const FileMapImageIteratorBase & r) : x(r.x), y(r.y), m_pImage(r.m_pImage)
	{
    }

    FileMapImageIteratorBase& operator=(const FileMapImageIteratorBase &r)
	{
		if(this!=&r)
		{
			x = r.x;
			y = r.y;
			m_pImage = r.m_pImage;
		}
        return *this;
    }
	image_type* m_pImage;
};

template <class PIXELTYPE> 
class FileMapImage;

template <class PIXELTYPE>
class FileMapImageIterator : public FileMapImageIteratorBase<FileMapImageIterator<PIXELTYPE>, FileMapImage<PIXELTYPE>, PIXELTYPE, PIXELTYPE &, PIXELTYPE *>
{
public:
	typedef FileMapImageIteratorBase<FileMapImageIterator, FileMapImage<PIXELTYPE>, PIXELTYPE, PIXELTYPE &, PIXELTYPE *> Base;

    FileMapImageIterator(int x = 0, int y = 0, FileMapImage<PIXELTYPE> * i = NULL)
    : Base(x, y, i)
    {}
};

template <class PIXELTYPE>
class ConstFileMapImageIterator : public FileMapImageIteratorBase<ConstFileMapImageIterator<PIXELTYPE>, const FileMapImage<PIXELTYPE>, PIXELTYPE, const PIXELTYPE &, const PIXELTYPE *>
{
public:
	typedef FileMapImageIteratorBase<ConstFileMapImageIterator, const FileMapImage<PIXELTYPE>, PIXELTYPE, const PIXELTYPE &, const PIXELTYPE *> Base;

    ConstFileMapImageIterator(int x = 0, int y = 0, const FileMapImage<PIXELTYPE> * i = NULL)
    : Base(x, y, i)
    {}
};

template <class Iterator>
class FileMapImageSequentialAccessIteratorPolicy
{
public:
    typedef Iterator BaseType;
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::difference_type::MoveX difference_type;
    typedef typename Iterator::reference reference;
    typedef typename Iterator::index_reference index_reference;
    typedef typename Iterator::pointer pointer;
	typedef typename Iterator::iterator_category iterator_category;

    static void initialize(BaseType & d) { }

    static reference dereference(const BaseType & d) {
        return *d;
    }

    static index_reference dereference(const BaseType & d, difference_type n) {
        int width = d.i->width();
        int dy = n / width;
        int dx = n % width;
        if (d.x + dx >= width) 
		{
			dy++; 
			dx -= width;
		}
        else if (d.x + dx < 0) 
		{
			dy--;
			dx += width;
		}
        return d(dx, dy);
    }

    static bool equal(const BaseType & d1, const BaseType & d2) {
        int width1 = d1.m_pImage->width();
        int width2 = d2.m_pImage->width();
        return (d1.y*width1 + d1.x) == (d2.y*width2 + d2.x);
    }

    static bool less(const BaseType & d1, const BaseType & d2) {
        int width1 = d1.m_pImage->width();
        int width2 = d2.m_pImage->width();
        return (d1.y*width1 + d1.x) < (d2.y*width2 + d2.x);
    }

    static difference_type difference(const BaseType & d1, const BaseType & d2) {
        int width1 = d1.m_pImage->width();
        int width2 = d2.m_pImage->width();
        return (d1.y*width1 + d1.x) - (d2.y*width2 + d2.x);
    }

    static void increment(BaseType & d) {
        ++d.x;
        if (d.x == d.m_pImage->width()) {
            d.x = 0;
            ++d.y;
        }
    }

    static void decrement(BaseType & d) {
        --d.x;
        if (d.x < 0) {
            d.x = d.m_pImage->width() - 1;
            --d.y;
        }
    }

    static void advance(BaseType & d, difference_type n) {
        int width = d.m_pImage->width();
        int dy = n / width;
        int dx = n % width;
        d.x += dx;
        d.y += dy;
        if (d.x >= width) 
		{
			++d.y; 
			d.x -= width;
		}
        if (d.x < 0) 
		{
			--d.y;
			d.x += width;
		}
    }

};

template <class PIXELTYPE>
class FileMapImage
{
public:
	typedef PIXELTYPE value_type;
    typedef PIXELTYPE PixelType;
	typedef PIXELTYPE & reference;
	typedef const PIXELTYPE & const_reference;
	typedef PIXELTYPE * pointer;
	typedef const PIXELTYPE * const_pointer;
	typedef FileMapImageIterator<PIXELTYPE> traverser;
	typedef ConstFileMapImageIterator<PIXELTYPE> const_traverser;
	typedef IteratorAdaptor<FileMapImageSequentialAccessIteratorPolicy<traverser> > iterator;
    typedef IteratorAdaptor<FileMapImageSequentialAccessIteratorPolicy<const_traverser> > const_iterator;
	typedef Diff2D difference_type;
	typedef Size2D size_type;
	typedef typename IteratorTraits<traverser>::DefaultAccessor Accessor;
	typedef typename IteratorTraits<const_traverser>::DefaultAccessor ConstAccessor;

	explicit FileMapImage(int w = 1, int h = 1, const PIXELTYPE & value = PIXELTYPE())
	{
		init(w,h,value);
	}

	explicit FileMapImage(const difference_type & size, const PIXELTYPE & value = PIXELTYPE())
	{
		init(size.x,size.y,value);
    }

	FileMapImage(const FileMapImage & image);
	~FileMapImage();

	FileMapImage & operator=(const FileMapImage & g);
	void resize(int w, int h, const PIXELTYPE & value = PIXELTYPE());

	int width() const
	{
		return m_width;
	}

	int height() const
	{
		return m_height;
	}

	size_type size() const
	{
		return size_type(width(), height());
	}

	reference operator[](difference_type const & d)
	{
		return getLine(d.y)[d.x];
	}

    const_reference operator[](difference_type const & d) const
	{
		return getLine(d.y)[d.x];
	}

	reference operator()(int dx, int dy)
	{
		return getLine(dy)[dx];
	}

    const_reference operator()(int dx, int dy) const
	{
		return getLine(dy)[dx];
	}

    pointer operator[](int dy)
	{
		return getLine(dy);
	}

    const_pointer operator[](int dy) const
	{
		return getLine(dy);
	}

	traverser upperLeft()
	{
		assert(m_width>0&&m_height>0);
		return traverser(0,0,this);
	}
    traverser lowerRight()
	{
		assert(m_width>0&&m_height>0);
		return traverser(m_width,m_height,this);
	}

    const_traverser upperLeft() const
	{
		assert(m_width>0&&m_height>0);
		return const_traverser(0,0,this);
	}
    const_traverser lowerRight() const
	{
		assert(m_width>0&&m_height>0);
		return const_traverser(m_width,m_height,this);
	}

    iterator begin()
	{
		assert(m_width>0&&m_height>0);
		return iterator(traverser(0,0,this));
	}
    iterator end()
	{
		assert(m_width>0&&m_height>0);
		return iterator(traverser(0,m_height,this));
	}

    const_iterator begin() const
	{
		assert(m_width>0&&m_height>0);
		return const_iterator(const_traverser(0,0,this));
	}

    const_iterator end() const
	{
		assert(m_width>0&&m_height>0);
		return const_iterator(const_traverser(0,m_height,this));
	}

	Accessor accessor() 
	{
        return Accessor();
    }

    ConstAccessor accessor() const 
	{
        return ConstAccessor();
    }

	void Destroy();

private:
	void switchBlock( int index ) const;
	void init(int w, int h, const PIXELTYPE & value = PIXELTYPE());
	PIXELTYPE* getLine(int line) const;
	void release();

	int m_width;
	int m_height;

#ifdef _WIN32
	HANDLE m_hFileMapping;
#else
	int m_fd;
#endif
	size_t m_block_size;
	size_t m_block_num;
	int m_lines_per_block;

	mutable PIXELTYPE ** m_lines;
	mutable int m_current_block;
	mutable unsigned char* m_pBlock;
};

template <class PIXELTYPE>
void FileMapImage<PIXELTYPE>::Destroy()
{
	//added by horace: To verify the correct destruction of FMImages.
	printf("\t\t release filemap image...\n");
	release();
}

template <class PIXELTYPE>
FileMapImage<PIXELTYPE>::FileMapImage(const FileMapImage & g)
{
	//这个函数的实现还可以优化，在init过程中把值设置好
	init(g.width(),g.height());
	const_iterator is = g.begin();
	const_iterator iend = g.end();
	iterator id = begin();
	for(;is!=iend;)
	{
		*(id++) = *(is++);
	}
}

template <class PIXELTYPE>
FileMapImage<PIXELTYPE>::~FileMapImage()
{
	//added by horace: To verify the correct destruction of FMImages.
	printf("\t\t release filemap image...\n");
	release();
}

template <class PIXELTYPE>
void FileMapImage<PIXELTYPE>::init(int w, int h, const PIXELTYPE & value)
{
	assert(w > 0 && h > 0);
	m_width = w;
	m_height = h;

#ifdef _WIN32
	SYSTEM_INFO sinf;
	GetSystemInfo(&sinf);
	size_t page_size = sinf.dwAllocationGranularity;
#else
	size_t page_size = getpagesize();
#endif

	size_t value_size = sizeof(value);
	int size_per_line = w * value_size;

	const int LESS = 3;

	m_block_size = (size_t)ceil( (double)(size_per_line<<LESS) / page_size ) * page_size;
	int lines_per_block = m_block_size/size_per_line;
	int lines_per_block_log2 = (int)floor( log((double)lines_per_block)/log(2.0) );
	m_lines_per_block = 1<<lines_per_block_log2;
	m_block_num = (size_t)ceil( (double)h/m_lines_per_block );

	char filenameTemplate[] = ".graph_tmpXXXXXX";

#ifdef _WIN32
	char *tmpReturn = mktemp(filenameTemplate);
	assert(tmpReturn);

	HANDLE hFile = CreateFileA(filenameTemplate,
		GENERIC_WRITE|GENERIC_READ,
		FILE_SHARE_READ|FILE_SHARE_WRITE,
		NULL,
		CREATE_ALWAYS,
		FILE_FLAG_DELETE_ON_CLOSE,
		NULL);

	assert(hFile!=INVALID_HANDLE_VALUE);

	//目前只处理最大int32能表示的最大值的文件内存映射，暂时不考虑使用高位，如果以后需要用到，这里的实现需要修改。
	m_hFileMapping = CreateFileMapping(hFile, NULL, PAGE_READWRITE, 0, m_block_size*m_block_num, NULL);
	CloseHandle(hFile);
	assert(m_hFileMapping);

	for(int i = 0; i < m_block_num; i++)
	{
		m_pBlock = (unsigned char*)MapViewOfFile(m_hFileMapping,FILE_SHARE_READ|FILE_SHARE_WRITE,0,i*m_block_size,m_block_size); 
		PIXELTYPE* pTmp = reinterpret_cast<PIXELTYPE*>(m_pBlock);
		std::uninitialized_fill_n(pTmp, m_lines_per_block*m_width, value);
		//FlushViewOfFile(m_pBlock, m_block_size);
		UnmapViewOfFile(m_pBlock);
	}

	m_pBlock  = (unsigned char*)MapViewOfFile(m_hFileMapping,FILE_SHARE_READ|FILE_SHARE_WRITE,0,0,m_block_size);
#else
	m_fd = mkstemp(filenameTemplate);
	assert(m_fd != -1);
	int nSuc = ftruncate(m_fd,m_block_size*m_block_num);
	assert(nSuc != -1);
	unlink(filenameTemplate);

	for(int i = 0; i < m_block_num; i++)
	{
		m_pBlock = (unsigned char*)mmap(NULL, m_block_size, PROT_READ|PROT_WRITE, MAP_SHARED, m_fd, i*m_block_size);
		PIXELTYPE* pTmp = reinterpret_cast<PIXELTYPE*>(m_pBlock);
		std::uninitialized_fill_n(pTmp, m_lines_per_block*m_width, value);
		munmap(m_pBlock, m_block_size);
	}

	m_pBlock = (unsigned char*)mmap(NULL, m_block_size, PROT_READ|PROT_WRITE, MAP_SHARED, m_fd, 0);
#endif

	m_current_block = 0;

	m_lines = new PIXELTYPE*[m_height];

	for (int i = 0; i < m_height; i++)
	{
		int block_index = i/m_lines_per_block;
		int line_offset = i%m_lines_per_block;
		if (block_index == m_current_block)
		{
			m_lines[i] = (PIXELTYPE*)(m_pBlock + line_offset*m_width*value_size);
		}
		else
		{
			m_lines[i] = NULL;
		}
	}
}

template <class PIXELTYPE>
PIXELTYPE* FileMapImage<PIXELTYPE>::getLine(int line) const
{
	assert(line>=0&&line<m_height);

	int block_index = line/m_lines_per_block;

	if (block_index != m_current_block)
	{
		switchBlock(block_index);
	}

	return m_lines[line];
}

template <class PIXELTYPE>
void FileMapImage<PIXELTYPE>::switchBlock(int index) const
{
	assert(index>=0&&index<m_block_num);
	assert(index!=m_current_block);

	int offset = m_current_block * m_lines_per_block;
	for (int i = 0; i < m_lines_per_block; i++)
	{
		int icurrent = offset + i;
		if (icurrent < m_height)
		{
			m_lines[icurrent] = NULL;
		}
	}
#ifdef _WIN32
	//FlushViewOfFile(m_pBlock, m_block_size);
	UnmapViewOfFile(m_pBlock);
	m_pBlock = (unsigned char*)MapViewOfFile(m_hFileMapping,FILE_SHARE_READ|FILE_SHARE_WRITE,0,m_block_size*index,m_block_size);
#else
	munmap(m_pBlock, m_block_size);
	m_pBlock = (unsigned char*)mmap(NULL, m_block_size, PROT_READ|PROT_WRITE, MAP_SHARED, m_fd, m_block_size*index);
#endif
	assert(m_pBlock);
	m_current_block = index;

	size_t size = sizeof(PIXELTYPE);

	offset = m_current_block * m_lines_per_block;
	for (int i=0; i < m_lines_per_block; i++)
	{
		int icurrent = offset + i;
		if (icurrent < m_height)
		{
			m_lines[icurrent] = (PIXELTYPE*)(m_pBlock + i*m_width*size);
		}
	}
}

template <class PIXELTYPE>
FileMapImage<PIXELTYPE> & FileMapImage<PIXELTYPE>::operator=(const FileMapImage<PIXELTYPE> & g)
{
	if (this!=&g)
	{
		if (this->width()!=g.width() || this->height()!=g.height())
		{
			resize(g.width(), g.height());
		}
		const_iterator is = g.begin();
		const_iterator iend = g.end();
		iterator id = begin();
		for(;is!=iend;)
		{
			*(id++) = *(is++);
		}
	}
	return *this;
}

template <class PIXELTYPE>
void FileMapImage<PIXELTYPE>::resize(int w, int h, const PIXELTYPE & value)
{
	assert(w>0&&h>0);
	release();
	init(w, h, value);
}

template <class PIXELTYPE>
void FileMapImage<PIXELTYPE>::release()
{
	if (m_lines)
	{
		delete[] m_lines;
		m_lines = NULL;
	}
#ifdef _WIN32
	//FlushViewOfFile(m_pBlock, m_block_size);
	UnmapViewOfFile(m_pBlock);
	CloseHandle(m_hFileMapping);
#else
	munmap(m_pBlock, m_block_size);
	close(m_fd);
#endif
}

template <class PIXELTYPE>
struct IteratorTraits<FileMapImageIterator<PIXELTYPE> >
{
    typedef FileMapImageIterator<PIXELTYPE>           Iterator;
    typedef Iterator                             iterator;
    typedef typename iterator::iterator_category iterator_category;
    typedef typename iterator::value_type        value_type;
    typedef typename iterator::reference         reference;
    typedef typename iterator::index_reference   index_reference;
    typedef typename iterator::pointer           pointer;
    typedef typename iterator::difference_type   difference_type;
    typedef typename iterator::row_iterator      row_iterator;
    typedef typename iterator::column_iterator   column_iterator;
	typedef typename AccessorTraits<PIXELTYPE>::default_accessor DefaultAccessor;
   // typedef StandardValueAccessor<PIXELTYPE>             DefaultAccessor;
    typedef DefaultAccessor             default_accessor;
};

template <class PIXELTYPE>
struct IteratorTraits<ConstFileMapImageIterator<PIXELTYPE> >
{
    typedef ConstFileMapImageIterator<PIXELTYPE>        Iterator;
    typedef Iterator                               iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
	typedef typename AccessorTraits<PIXELTYPE>::default_const_accessor DefaultAccessor;
    //typedef StandardConstValueAccessor<PIXELTYPE>          DefaultAccessor;
    typedef DefaultAccessor          default_accessor;
};

/*
template <class PIXELTYPE>
struct IteratorTraits<FileMapImageIterator<RGBValue<PIXELTYPE> > >
{
    typedef FileMapImageIterator<RGBValue<PIXELTYPE> >           Iterator;
    typedef Iterator                             iterator;
    typedef typename iterator::iterator_category iterator_category;
    typedef typename iterator::value_type        value_type;
    typedef typename iterator::reference         reference;
    typedef typename iterator::index_reference   index_reference;
    typedef typename iterator::pointer           pointer;
    typedef typename iterator::difference_type   difference_type;
    typedef typename iterator::row_iterator      row_iterator;
    typedef typename iterator::column_iterator   column_iterator;
    typedef RGBAccessor<RGBValue<PIXELTYPE> >            DefaultAccessor;
    typedef RGBAccessor<RGBValue<PIXELTYPE> >            default_accessor;
};

template <class PIXELTYPE>
struct IteratorTraits<ConstFileMapImageIterator<RGBValue<PIXELTYPE> > >
{
    typedef ConstFileMapImageIterator<RGBValue<PIXELTYPE> >        Iterator;
    typedef Iterator                               iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
    typedef RGBAccessor<RGBValue<PIXELTYPE> >              DefaultAccessor;
    typedef RGBAccessor<RGBValue<PIXELTYPE> >              default_accessor;
};
*/

typedef FileMapImage<UInt8> UInt8FMImage;
typedef FileMapImage<Int8> Int8FMImage;
typedef FileMapImage<UInt16> UInt16FMImage;
typedef FileMapImage<Int16> Int16FMImage;
typedef FileMapImage<UInt32> UInt32FMImage;
typedef FileMapImage<Int32> Int32FMImage;
typedef FileMapImage<float> FFMImage;
typedef FileMapImage<double> DFMImage;
typedef FileMapImage<RGBValue<UInt8> > UInt8RGBFMImage;
typedef FileMapImage<RGBValue<Int8> > Int8RGBFMImage;
typedef FileMapImage<RGBValue<UInt16> > UInt16RGBFMImage;
typedef FileMapImage<RGBValue<Int16> > Int16RGBFMImage;
typedef FileMapImage<RGBValue<UInt32> > UInt32RGBFMImage;
typedef FileMapImage<RGBValue<Int32> > Int32RGBFMImage;
typedef FileMapImage<RGBValue<float> > FRGBFMImage;
typedef FileMapImage<RGBValue<double> > DRGBFMImage;


template <class PixelType>
inline triple<
	typename FileMapImage<PixelType>::const_traverser,
	typename FileMapImage<PixelType>::const_traverser,
	typename FileMapImage<PixelType>::ConstAccessor
>
srcImageRange(FileMapImage<PixelType> const & img, Rect2D const & roi)
{
	return triple<
		typename FileMapImage<PixelType>::const_traverser,
		typename FileMapImage<PixelType>::const_traverser, 
		typename FileMapImage<PixelType>::ConstAccessor>
		(img.upperLeft()+roi.upperLeft(),
		img.upperLeft() + roi.lowerRight(),
		img.accessor());
}

template <class PixelType, class Accessor>
inline triple<typename FileMapImage<PixelType>::const_traverser, 
              typename FileMapImage<PixelType>::const_traverser, Accessor>
srcImageRange(FileMapImage<PixelType> const & img, Accessor a)
{
    return triple<typename FileMapImage<PixelType>::const_traverser, 
                  typename FileMapImage<PixelType>::const_traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename FileMapImage<PixelType>::const_traverser, Accessor>
srcImage(FileMapImage<PixelType> const & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}


template <class PixelType>
inline triple<
    typename FileMapImage<PixelType>::traverser, 
	typename FileMapImage<PixelType>::traverser,
	typename FileMapImage<PixelType>::Accessor
>
destImageRange(FileMapImage<PixelType> & img, Rect2D const & roi)
{
	return triple<
		typename FileMapImage<PixelType>::traverser,
		typename FileMapImage<PixelType>::traverser, 
		typename FileMapImage<PixelType>::Accessor>
		(img.upperLeft()+roi.upperLeft(),
		img.upperLeft() + roi.lowerRight(),
		img.accessor());
}


template <class PixelType, class Accessor>
inline triple<typename FileMapImage<PixelType>::traverser, 
              typename FileMapImage<PixelType>::traverser, Accessor>
destImageRange(FileMapImage<PixelType> & img, Accessor a)
{
    return triple<typename FileMapImage<PixelType>::traverser, 
                  typename FileMapImage<PixelType>::traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename FileMapImage<PixelType>::traverser, Accessor>
destImage(FileMapImage<PixelType> & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType>::traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename FileMapImage<PixelType>::const_traverser, Accessor>
maskImage(FileMapImage<PixelType> const & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType>
inline triple<typename FileMapImage<PixelType>::const_traverser, 
              typename FileMapImage<PixelType>::const_traverser, 
              typename FileMapImage<PixelType>::ConstAccessor>
srcImageRange(FileMapImage<PixelType> const & img)
{
    return triple<typename FileMapImage<PixelType>::const_traverser, 
                  typename FileMapImage<PixelType>::const_traverser, 
                  typename FileMapImage<PixelType>::ConstAccessor>(img.upperLeft(),
                                                                 img.lowerRight(),
                                                                 img.accessor());
}

template <class PixelType>
inline pair< typename FileMapImage<PixelType>::const_traverser, 
             typename FileMapImage<PixelType>::ConstAccessor>
srcImage(FileMapImage<PixelType> const & img)
{
    return pair<typename FileMapImage<PixelType>::const_traverser, 
                typename FileMapImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

template <class PixelType>
inline triple< typename FileMapImage<PixelType>::traverser, 
               typename FileMapImage<PixelType>::traverser, 
               typename FileMapImage<PixelType>::Accessor>
destImageRange(FileMapImage<PixelType> & img)
{
    return triple<typename FileMapImage<PixelType>::traverser, 
                  typename FileMapImage<PixelType>::traverser, 
                  typename FileMapImage<PixelType>::Accessor>(img.upperLeft(),
                                                            img.lowerRight(),
                                                            img.accessor());
}

template <class PixelType>
inline pair< typename FileMapImage<PixelType>::traverser, 
             typename FileMapImage<PixelType>::Accessor>
destImage(FileMapImage<PixelType> & img)
{
    return pair<typename FileMapImage<PixelType>::traverser, 
                typename FileMapImage<PixelType>::Accessor>(img.upperLeft(), 
                                                          img.accessor());
}

template <class PixelType>
inline pair< typename FileMapImage<PixelType>::const_traverser, 
             typename FileMapImage<PixelType>::ConstAccessor>
maskImage(FileMapImage<PixelType> const & img)
{
    return pair<typename FileMapImage<PixelType>::const_traverser, 
                typename FileMapImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

} // namespace vigra

using namespace vigra;

#endif