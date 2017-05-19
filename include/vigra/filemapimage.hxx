#ifndef DEFIER_GRAPH_H
#define DEFIER_GRAPH_H

#include <io.h>

#include <assert.h>

#ifndef _WIN32
#include <unistd.h>
#include <sys/mman.h>
#else
#include <windows.h>
#undef max
#undef min
#endif

#pragma push_macro("DIFFERENCE")
#undef DIFFERENCE
#include <vigra/diff2d.hxx>
#include <vigra/rgbvalue.hxx>
#pragma pop_macro("DIFFERENCE")

#include <list>

using std::list;

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

template <class PIXELTYPE, class IMAGETYPE>
class FileMapImageIterator : public FileMapImageIteratorBase<FileMapImageIterator<PIXELTYPE, IMAGETYPE>, IMAGETYPE, PIXELTYPE, PIXELTYPE &, PIXELTYPE *>
{
public:
	typedef FileMapImageIteratorBase<FileMapImageIterator, IMAGETYPE, PIXELTYPE, PIXELTYPE &, PIXELTYPE *> Base;

    FileMapImageIterator(int x = 0, int y = 0, IMAGETYPE * i = NULL)
    : Base(x, y, i)
    {}
};

template <class PIXELTYPE, class IMAGETYPE>
class ConstFileMapImageIterator : public FileMapImageIteratorBase<ConstFileMapImageIterator<PIXELTYPE, IMAGETYPE>, const IMAGETYPE, PIXELTYPE, const PIXELTYPE &, const PIXELTYPE *>
{
public:
	typedef FileMapImageIteratorBase<ConstFileMapImageIterator, const IMAGETYPE, PIXELTYPE, const PIXELTYPE &, const PIXELTYPE *> Base;

    ConstFileMapImageIterator(int x = 0, int y = 0, const IMAGETYPE * i = NULL)
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

template <class PIXELTYPE, class Alloc = std::allocator<PIXELTYPE> >
class FileMapImage
{
public:
	typedef FileMapImage<PIXELTYPE, Alloc> self_type;

	typedef PIXELTYPE value_type;
    typedef PIXELTYPE PixelType;
	typedef PIXELTYPE & reference;
	typedef const PIXELTYPE & const_reference;
	typedef PIXELTYPE * pointer;
	typedef const PIXELTYPE * const_pointer;

	typedef FileMapImageIterator<PIXELTYPE, self_type> traverser;
	typedef ConstFileMapImageIterator<PIXELTYPE, self_type> const_traverser;

	typedef FileMapImageIterator<PIXELTYPE, self_type> Iterator;
	typedef ConstFileMapImageIterator<PIXELTYPE, self_type> ConstIterator;

	typedef IteratorAdaptor<FileMapImageSequentialAccessIteratorPolicy<traverser> > iterator;
    typedef IteratorAdaptor<FileMapImageSequentialAccessIteratorPolicy<const_traverser> > const_iterator;
	typedef Diff2D difference_type;
	typedef Size2D size_type;
	typedef typename IteratorTraits<traverser>::DefaultAccessor Accessor;
	typedef typename IteratorTraits<const_traverser>::DefaultAccessor ConstAccessor;

	typedef Alloc allocator_type;

	typedef Alloc Allocator;
	typedef typename Alloc::template rebind<PIXELTYPE *>::other LineAllocator;

	FileMapImage(int w = 1, int h = 1, const value_type & value = value_type(), const Alloc & alloc = Alloc())
		:m_allocator(alloc)
		,m_pallocator(alloc)
	{
		initMembers();
		resize(w, h, value);
	}

	explicit FileMapImage(const difference_type & size, const value_type & value = value_type(), const Alloc & alloc = Alloc())
		:m_allocator(alloc)
		,m_pallocator(alloc)
	{
		initMembers();
		resize(size.x, size.y, value);
    }

	FileMapImage(const FileMapImage & rhs)
	{
		initMembers();
		resize(rhs.width(), rhs.height(), value_type());

		const_iterator is = rhs.begin();
        const_iterator iend = rhs.end();
        iterator id = begin();

		while(is != iend)
		{
			*(id++) = *(is++);
		}
	}

	~FileMapImage()
	{
		closeTmpFile();
		deallocate();

		if(m_blockIsClean)
		{
			delete[] m_blockIsClean;
		}

		if(m_blockInFile)
		{
			delete[] m_blockInFile;
		}

		if(m_blocksInMemory)
		{
			delete m_blocksInMemory;
		}
	}

	FileMapImage & operator=(const FileMapImage & rhs)
	{
		if(this != &rhs)
		{
			resize(rhs.width(), rhs.height());

			const_iterator is = rhs.begin();
			const_iterator iend = rhs.end();
			iterator id = begin();

			while(is != iend)
			{
				*(id++) = *(is++);
			}
		}
		return *this;
	}

	FileMapImage & operator=(value_type pixel)
	{
		iterator i = begin();
		iterator iend = end();

		while(i != iend)
		{
			*(i++) = pixel;
		}
		return *this;
	}

	void resize(int width, int height)
	{
		if(width != m_width || height != m_height)
		{
            resize(width, height, value_type());
		}
	}

	void resize(const difference_type & size)
	{
		if(size.x != m_width || size.y != m_height)
        {
            resize(size.x, size.y, value_type());
        }
	}

	void resize(int width, int height, const value_type & d)
	{
		resizeImpl(width, height, d);
	}

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

	bool isInside(const difference_type & d) const
	{
		return d.x >= 0 && d.y >=0 && d.x < m_width && d.y < m_height;
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

private:

	void deallocate();
	void closeTmpFile();
	void resizeImpl(int width, int height, const value_type & d);

	PIXELTYPE* getLine(int y) const;
	PIXELTYPE* getLine(int y);
	PIXELTYPE* getLineMiss(int y) const;

	PIXELTYPE* swapOutBlock(int index) const;
	void loadInBlock(int index, PIXELTYPE* pBlock) const;

	void initTmpfile() const;

	inline int block_max_size()
	{
		return 1<<20;
	}

	inline int block_max_num()
	{
		return 32;
	}

	void initMembers()
	{
		m_blocksNeeded = 0;
		m_blocksMaxLoadNum = 0;
		m_linesPerBlock = 0;
		m_linesPerBlockLog2 = 0;

		m_width = 0;
		m_height = 0;
		m_initPixel = value_type();

		m_data = 0;
		m_lines = 0;

		m_blockIsClean = 0;
		m_blockInFile = 0;
		m_blocksInMemory = 0;

		m_tmpFilename = 0;

#ifdef _WIN32
		m_hFile = INVALID_HANDLE_VALUE;
#else
		m_file = 0;
#endif
	}

	int m_blocksNeeded;
	int m_blocksMaxLoadNum;
	int m_linesPerBlock;
	int m_linesPerBlockLog2;
	
	int m_width;
	int m_height;
	PIXELTYPE m_initPixel;

	Alloc m_allocator;
	LineAllocator m_pallocator;

	mutable PIXELTYPE * m_data;
	mutable PIXELTYPE ** m_lines;
	mutable bool * m_blockIsClean;
	mutable bool * m_blockInFile;
	mutable list<int>* m_blocksInMemory;
	mutable list<int>::iterator m_mostRecentlyLoadedBlockIterator;

	mutable char* m_tmpFilename;
#ifdef _WIN32
	mutable HANDLE m_hFile;
#else
	mutable FILE * m_file;
#endif
};

template<class PIXELTYPE, class Alloc>
void FileMapImage<PIXELTYPE, Alloc>::resizeImpl(int width, int height, const value_type & d)
{
	assert(width > 0 && height > 0);
	
	const int BLOCK_MAX_SIZE = block_max_size();
	const int BLOCK_MAX_NUM = block_max_num();

	int size_of_pixsel = sizeof(value_type);

	int linesPerBlock = (int)floor( (double)BLOCK_MAX_SIZE / (width * size_of_pixsel) );

	assert(linesPerBlock > 0);

	int linesPerBlockLog2 = (int)floor( log((double)linesPerBlock) / log(2.0) );
	linesPerBlock = 1<<linesPerBlockLog2;
	int blocksNeeded = (int)ceil( (double)height / linesPerBlock );
	int blocksMaxLoadNum = blocksNeeded < BLOCK_MAX_NUM ? blocksNeeded : BLOCK_MAX_NUM;

	int old_size = m_blocksMaxLoadNum * m_linesPerBlock * m_width;
	int new_size = blocksMaxLoadNum * linesPerBlock * width;

	if(width != m_width || height != m_height)
	{
		value_type * newdata = 0;
		value_type ** newlines = 0;

		if(new_size!=old_size)
		{
			newdata = m_allocator.allocate( typename Alloc::size_type(new_size) );
			std::uninitialized_fill_n(newdata, new_size, d);
			newlines = m_pallocator.allocate( typename Alloc::size_type(height) );
			std::uninitialized_fill_n(newlines, height, (value_type*)0);
			deallocate();
		}
		else
		{
			newdata = m_data;
			std::fill_n(newdata, new_size, d);
			newlines = m_pallocator.allocate( typename Alloc::size_type(height) );
			std::uninitialized_fill_n(newlines, height, (value_type*)0);
			m_pallocator.deallocate(m_lines, typename Alloc::size_type(m_height));
		}

		m_data = newdata;
        m_lines = newlines;
        m_width = width;
        m_height = height;
	}
	else
	{
		std::fill_n(m_data, new_size, d);
		std::fill_n(m_lines, m_height, (value_type*)0);
	}
	
	m_initPixel = d;

	m_blocksNeeded = blocksNeeded;
	m_blocksMaxLoadNum = blocksMaxLoadNum;
	m_linesPerBlock = linesPerBlock;
	m_linesPerBlockLog2 = linesPerBlockLog2;

	if(m_blockIsClean)
	{
		delete[] m_blockIsClean;
		m_blockIsClean = 0;
	}

	if(m_blockInFile)
	{
		delete[] m_blockInFile;
		m_blockInFile = 0;
	}

	m_blockIsClean = new bool[m_blocksNeeded];
	m_blockInFile = new bool[m_blocksNeeded];

	std::fill_n(m_blockIsClean, m_blocksNeeded, true);
	std::fill_n(m_blockInFile, m_blocksNeeded, false);

	if(!m_blocksInMemory)
	{
		m_blocksInMemory = new list<int>();
	}
	m_blocksInMemory->clear();

	m_mostRecentlyLoadedBlockIterator = m_blocksInMemory->begin();

	closeTmpFile();
}

template<class PIXELTYPE, class Alloc>
void FileMapImage<PIXELTYPE, Alloc>::deallocate()
{
	if(m_data)
    {
		for(int i = 0; i < m_blocksMaxLoadNum * m_linesPerBlock * m_width; i++)
		{
			m_data[i].~PIXELTYPE();
		}

        m_allocator.deallocate(m_data, typename Alloc::size_type(m_blocksMaxLoadNum * m_linesPerBlock * m_width));
        m_pallocator.deallocate(m_lines, typename Alloc::size_type(m_height));
    }
}

template<class PIXELTYPE, class Alloc>
void FileMapImage<PIXELTYPE, Alloc>::closeTmpFile()
{
#ifdef _WIN32
	if(m_hFile != INVALID_HANDLE_VALUE)
	{
		CloseHandle(m_hFile);
		m_hFile = INVALID_HANDLE_VALUE;
	}
#else
	if(m_file)
	{
		fclose(m_file);
		m_file = 0;
	}
#endif

	if(m_tmpFilename)
	{
		delete[] m_tmpFilename;
		m_tmpFilename = 0;
	}
}

template <class PIXELTYPE, class Alloc>
PIXELTYPE* FileMapImage<PIXELTYPE, Alloc>::getLine(int y) const
{
	PIXELTYPE* line = m_lines[y];

	if(!line)
	{
		line = getLineMiss(y);
	}

	return line;
}

template <class PIXELTYPE, class Alloc>
PIXELTYPE* FileMapImage<PIXELTYPE, Alloc>::getLine(int y)
{
	PIXELTYPE* line = m_lines[y];

	if(!line)
	{
		line = getLineMiss(y);
	}

	m_blockIsClean[y>>m_linesPerBlockLog2] = false;

	return line;
}

template <class PIXELTYPE, class Alloc>
PIXELTYPE* FileMapImage<PIXELTYPE, Alloc>::getLineMiss(int y) const
{
	int blockNum = y>>m_linesPerBlockLog2;

	int from = blockNum * m_linesPerBlock;
	int to = from + m_linesPerBlock;

	if(to > m_height)
	{
		to = m_height;
	}

	PIXELTYPE* pBlock = 0;

	int size = m_blocksInMemory->size();

	if(size < m_blocksMaxLoadNum)
	{
		pBlock = m_data + m_linesPerBlock * size  * m_width;
	}
	else
	{
		int blockNumForSwapOut = 0;
		if(m_mostRecentlyLoadedBlockIterator == m_blocksInMemory->begin())
		{
			blockNumForSwapOut = m_blocksInMemory->back();
            m_blocksInMemory->pop_back();
		}
		else
		{
			list<int>::iterator candidate = m_mostRecentlyLoadedBlockIterator;
			--candidate;
			blockNumForSwapOut = *candidate;
			m_blocksInMemory->erase(candidate);
		}
		pBlock = swapOutBlock(blockNumForSwapOut);
	}

	loadInBlock(blockNum, pBlock);

	m_blockIsClean[blockNum] = true;

    list<int>::iterator i = m_blocksInMemory->begin();
    for(; i != m_blocksInMemory->end(); i++)
	{
		if(*i > blockNum)
		{
			break;
		}
    }
    m_blocksInMemory->insert(i, blockNum);

    m_mostRecentlyLoadedBlockIterator = --i;
    for(i++; i != m_blocksInMemory->end(); i++)
	{
		if(*i == ++blockNum)
		{
			m_mostRecentlyLoadedBlockIterator = i;
		}
        else
		{
			break;
		}
    }

	return m_lines[y];
}

template <class PIXELTYPE, class Alloc>
PIXELTYPE* FileMapImage<PIXELTYPE, Alloc>::swapOutBlock(int index) const
{
	int from = index * m_linesPerBlock;
	int to = from + m_linesPerBlock;

	if(to > m_height)
	{
		to = m_height;
	}

	PIXELTYPE* pBlock = m_lines[from];

	std::fill_n(m_lines + from, to - from, (value_type*)0);

	if(!m_blockIsClean[index])
	{
		m_blockInFile[index] = true;
		int pixelsToWrite = (to - from) * m_width;

#ifdef _WIN32
		if(m_hFile == INVALID_HANDLE_VALUE)
#else
		if(!m_file) 
#endif
		{
			initTmpfile();
		}

#ifdef _WIN32
		DWORD dwError = NO_ERROR;
		LONGLONG offset = (LONGLONG)m_width
			* (LONGLONG)from
			* (LONGLONG)sizeof(PIXELTYPE);
		LARGE_INTEGER liOffset;
		liOffset.QuadPart = offset;
		liOffset.LowPart = SetFilePointer(m_hFile, liOffset.LowPart, &liOffset.HighPart, FILE_BEGIN);
		if(liOffset.LowPart == INVALID_SET_FILE_POINTER && (dwError = GetLastError()) != NO_ERROR)
		{
			LPVOID lpMsgBuf;
			FormatMessage(
				FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
				NULL,
				dwError,
				MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
				(LPTSTR) &lpMsgBuf,
				0,
				NULL);
			cerr << endl << lpMsgBuf << endl;
			LocalFree(lpMsgBuf);
			vigra_fail("Unable to seek within temporary file.\n");
		}

		DWORD bytesWritten;
		if(0 == WriteFile(m_hFile, pBlock, sizeof(PIXELTYPE) * pixelsToWrite, &bytesWritten, NULL))
		{
			DWORD dwError = GetLastError();
			LPVOID lpMsgBuf;
			FormatMessage(
				FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
				NULL,
				dwError,
				MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
				(LPTSTR) &lpMsgBuf,
				0,
				NULL);
			cerr << endl << lpMsgBuf << endl;
			LocalFree(lpMsgBuf);
			vigra_fail("pwblend: error writing to image swap file.\n");
		}
#else
		off_t offset = (off_t)m_width
			* (off_t)from
			* (off_t)sizeof(PIXELTYPE);
		if(fseeko(m_file, offset, SEEK_SET) != 0)
		{
			vigra_fail(strerror(errno));
		}

		clearerr(m_file);
		int itemsWritten = fwrite(pBlock, sizeof(PIXELTYPE), pixelsToWrite, m_file);
		if (itemsWritten < pixelsToWrite)
		{
			perror("pwblend");
			vigra_fail("pwblend: error writing to image swap file.\n");
		}
#endif
	}

	return pBlock;
}

template <class PIXELTYPE, class Alloc>
void FileMapImage<PIXELTYPE, Alloc>::loadInBlock(int index, PIXELTYPE* pBlock) const
{
	int from = index * m_linesPerBlock;
	int to = from + m_linesPerBlock;

	if(to > m_height)
	{
		to = m_height;
	}

	for(int i = from; i < to; i++)
	{
		m_lines[i] = pBlock + (i - from) * m_width; 
	}

	int pixelsToRead = (to - from) * m_width;

	if(m_blockInFile[index])
	{
#ifdef _WIN32
		DWORD dwError = NO_ERROR;
		LONGLONG offset = (LONGLONG)m_width
			* (LONGLONG)from
			* (LONGLONG)sizeof(PIXELTYPE);
		LARGE_INTEGER liOffset;
		liOffset.QuadPart = offset;
		liOffset.LowPart = SetFilePointer(m_hFile, liOffset.LowPart, &liOffset.HighPart, FILE_BEGIN);
		if(liOffset.LowPart == INVALID_SET_FILE_POINTER && (dwError = GetLastError()) != NO_ERROR)
		{
			LPVOID lpMsgBuf;
			FormatMessage(
				FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
				NULL,
				dwError,
				MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
				(LPTSTR) &lpMsgBuf,
				0,
				NULL);
			cerr << endl << lpMsgBuf << endl;
			LocalFree(lpMsgBuf);
			vigra_fail("Unable to seek within temporary file.\n");
		}

		DWORD bytesRead;
		if(0 == ReadFile(m_hFile, pBlock, sizeof(PIXELTYPE) * pixelsToRead, &bytesRead, NULL))
		{
			DWORD dwError = GetLastError();
			LPVOID lpMsgBuf;
			FormatMessage(
				FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
				NULL,
				dwError,
				MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
				(LPTSTR) &lpMsgBuf,
				0,
				NULL);
			cerr << endl << lpMsgBuf << endl;
			LocalFree(lpMsgBuf);
			vigra_fail("pwblend: error reading from image swap file.\n");
		}
#else
		off_t offset = (off_t)m_width
			* (off_t)from
			* (off_t)sizeof(PIXELTYPE); 
		if(fseeko(tmpFile_, offset, SEEK_SET) != 0)
		{
			vigra_fail(strerror(errno));
		}

		clearerr(m_file);
		int itemsRead = fread(pBlock, sizeof(PIXELTYPE), pixelsToRead, m_file);
		if (itemsRead < pixelsToRead) {
			perror("pwblend");
			vigra_fail("pwblend: error reading from image swap file.\n");
		}
#endif
	}
	else
	{
		std::fill_n(m_lines[from], pixelsToRead, m_initPixel);
	}
}

template <class PIXELTYPE, class Alloc>
void FileMapImage<PIXELTYPE, Alloc>::initTmpfile() const
{
	char filenameTemplate[] = ".graph_tmp00XXXXXX";
#ifdef _WIN32
	char *tmpReturn = mktemp(filenameTemplate);
	int n = 1;
	while(!tmpReturn && n < 100)
	{
		int a = n % 10;
		int b = (n-a) / 10;
		char tempbuf[] = ".graph_tmp00XXXXXX";
		tempbuf[10] += b;
		tempbuf[11] += a;
		tmpReturn = mktemp(tempbuf);
		if(tmpReturn)
		{
			strcpy(filenameTemplate, tempbuf);
			break;
		}
		n++;
	}
	assert(tmpReturn);

	m_hFile = CreateFileA(filenameTemplate,
		GENERIC_WRITE|GENERIC_READ,
		FILE_SHARE_READ|FILE_SHARE_WRITE,
		NULL,
		CREATE_ALWAYS,
		FILE_FLAG_DELETE_ON_CLOSE,
		NULL);

	assert(m_hFile!=INVALID_HANDLE_VALUE);
#else
	int tmpFD = mkstemp(filenameTemplate);
	int n = 1;
	while(tmpFD == -1 && n < 100)
	{
		int a = n % 10;
		int b = (n-a) / 10;
		char tempbuf[] = ".graph_tmp00XXXXXX";
		tempbuf[10] += b;
		tempbuf[11] += a;
		tmpFD = mkstemp(tempbuf);
		if(tempbuf != -1)
		{
			strcpy(filenameTemplate, tempbuf);
			break;
		}
	}
	assert(tmfFD != -1);

	m_file = fdopen(tmpFD, "wb+");
	assert(m_file);

	unlink(filenameTemplate);
#endif

	unsigned int filenameTemplateLength = (unsigned int)strlen(filenameTemplate) + 1;
	m_tmpFilename = new char[filenameTemplateLength];
	strncpy(m_tmpFilename, filenameTemplate, filenameTemplateLength);
}

template <class PIXELTYPE, class IMAGETYPE>
struct IteratorTraits<FileMapImageIterator<PIXELTYPE, IMAGETYPE> >
{
    typedef FileMapImageIterator<PIXELTYPE, IMAGETYPE>           Iterator;
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

template <class PIXELTYPE, class IMAGETYPE>
struct IteratorTraits<ConstFileMapImageIterator<PIXELTYPE, IMAGETYPE> >
{
    typedef ConstFileMapImageIterator<PIXELTYPE, IMAGETYPE>        Iterator;
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

template <class PixelType, class Accessor, class Alloc>
inline triple<
	typename FileMapImage<PixelType, Alloc>::const_traverser, 
	typename FileMapImage<PixelType, Alloc>::const_traverser,
	Accessor
>
srcImageRange(FileMapImage<PixelType, Alloc> const & img, Accessor a)
{
    return triple<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                  typename FileMapImage<PixelType, Alloc>::const_traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor, class Alloc>
inline triple<
	typename FileMapImage<PixelType, Alloc>::const_traverser,
	typename FileMapImage<PixelType, Alloc>::const_traverser,
	Accessor
>
srcImageRange(FileMapImage<PixelType, Alloc> const & img, Rect2D const & roi, Accessor a)
{
	vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
		roi.right() <= img.width() && roi.bottom() <= img.height(),
		"srcImageRange(): ROI rectangle outside image.");
	return triple<
		typename FileMapImage<PixelType, Alloc>::const_traverser,
		typename FileMapImage<PixelType, Alloc>::const_traverser, 
		Accessor
	>
	(
	img.upperLeft() + roi.upperLeft(),
	img.upperLeft() + roi.lowerRight(),
	a
	);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::const_traverser, Accessor>
srcImage(FileMapImage<PixelType, Alloc> const & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::const_traverser, Accessor>
srcImage(FileMapImage<PixelType, Alloc> const & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "srcImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser,
                Accessor>(img.upperLeft() + ul, a);
}

template <class PixelType, class Accessor, class Alloc>
inline triple<typename FileMapImage<PixelType>::traverser, 
              typename FileMapImage<PixelType>::traverser, Accessor>
destImageRange(FileMapImage<PixelType, Alloc> & img, Accessor a)
{
    return triple<typename FileMapImage<PixelType, Alloc>::traverser, 
                  typename FileMapImage<PixelType, Alloc>::traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor, class Alloc>
inline triple<
    typename FileMapImage<PixelType, Alloc>::traverser, 
	typename FileMapImage<PixelType, Alloc>::traverser,
	Accessor
>
destImageRange(FileMapImage<PixelType, Alloc> & img, Rect2D const & roi, Accessor a)
{
	vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "destImageRange(): ROI rectangle outside image.");
	return triple<
		typename FileMapImage<PixelType, Alloc>::traverser,
		typename FileMapImage<PixelType, Alloc>::traverser, 
		Accessor>
		(img.upperLeft() + roi.upperLeft(),
		img.upperLeft() + roi.lowerRight(),
		a);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::traverser, Accessor>
destImage(FileMapImage<PixelType, Alloc> & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType, Alloc>::traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::traverser, Accessor>
destImage(FileMapImage<PixelType, Alloc> & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "destImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::traverser,
                Accessor>(img.upperLeft() + ul, a);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::const_traverser, Accessor>
maskImage(FileMapImage<PixelType, Alloc> const & img, Accessor a)
{
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor, class Alloc>
inline pair<typename FileMapImage<PixelType, Alloc>::const_traverser, Accessor>
maskImage(FileMapImage<PixelType, Alloc> const & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "maskImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser,
                Accessor>(img.upperLeft() + ul, a);
}

template <class PixelType, class Alloc>
inline triple<
	typename FileMapImage<PixelType, Alloc>::const_traverser, 
	typename FileMapImage<PixelType, Alloc>::const_traverser, 
	typename FileMapImage<PixelType, Alloc>::ConstAccessor
>
srcImageRange(FileMapImage<PixelType, Alloc> const & img)
{
    return triple<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                  typename FileMapImage<PixelType, Alloc>::const_traverser, 
                  typename FileMapImage<PixelType, Alloc>::ConstAccessor>(img.upperLeft(),
                                                                 img.lowerRight(),
                                                                 img.accessor());
}

template <class PixelType, class Alloc>
inline triple<
	typename FileMapImage<PixelType, Alloc>::const_traverser,
	typename FileMapImage<PixelType, Alloc>::const_traverser,
	typename FileMapImage<PixelType, Alloc>::ConstAccessor
>
srcImageRange(FileMapImage<PixelType, Alloc> const & img, Rect2D const & roi)
{
	vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "srcImageRange(): ROI rectangle outside image.");
	return triple<
		typename FileMapImage<PixelType, Alloc>::const_traverser,
		typename FileMapImage<PixelType, Alloc>::const_traverser, 
		typename FileMapImage<PixelType, Alloc>::ConstAccessor
	>
	(
	img.upperLeft() + roi.upperLeft(),
	img.upperLeft() + roi.lowerRight(),
	img.accessor()
	);
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::const_traverser, 
             typename FileMapImage<PixelType, Alloc>::ConstAccessor>
srcImage(FileMapImage<PixelType, Alloc> const & img)
{
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                typename FileMapImage<PixelType, Alloc>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::const_traverser,
             typename FileMapImage<PixelType, Alloc>::ConstAccessor>
srcImage(FileMapImage<PixelType, Alloc> const & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "srcImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser,
                typename FileMapImage<PixelType, Alloc>::ConstAccessor>(img.upperLeft() + ul,
                                                                      img.accessor());
}

template <class PixelType, class Alloc>
inline triple< typename FileMapImage<PixelType, Alloc>::traverser, 
               typename FileMapImage<PixelType, Alloc>::traverser, 
               typename FileMapImage<PixelType, Alloc>::Accessor>
destImageRange(FileMapImage<PixelType, Alloc> & img)
{
    return triple<typename FileMapImage<PixelType, Alloc>::traverser, 
                  typename FileMapImage<PixelType, Alloc>::traverser, 
                  typename FileMapImage<PixelType, Alloc>::Accessor>(img.upperLeft(),
                                                            img.lowerRight(),
                                                            img.accessor());
}

template <class PixelType, class Alloc>
inline triple<
    typename FileMapImage<PixelType, Alloc>::traverser, 
	typename FileMapImage<PixelType, Alloc>::traverser,
	typename FileMapImage<PixelType, Alloc>::Accessor
>
destImageRange(FileMapImage<PixelType, Alloc> & img, Rect2D const & roi)
{
	vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "destImageRange(): ROI rectangle outside image.");
	return triple<
		typename FileMapImage<PixelType, Alloc>::traverser,
		typename FileMapImage<PixelType, Alloc>::traverser, 
		typename FileMapImage<PixelType, Alloc>::Accessor>
		(img.upperLeft() + roi.upperLeft(),
		img.upperLeft() + roi.lowerRight(),
		img.accessor());
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::traverser, 
             typename FileMapImage<PixelType, Alloc>::Accessor>
destImage(FileMapImage<PixelType, Alloc> & img)
{
    return pair<typename FileMapImage<PixelType, Alloc>::traverser, 
                typename FileMapImage<PixelType, Alloc>::Accessor>(img.upperLeft(), 
                                                          img.accessor());
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::traverser,
             typename FileMapImage<PixelType, Alloc>::Accessor>
destImage(FileMapImage<PixelType, Alloc> & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "destImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::traverser,
                typename FileMapImage<PixelType, Alloc>::Accessor>(img.upperLeft() + ul,
                                                                 img.accessor());
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::const_traverser, 
             typename FileMapImage<PixelType, Alloc>::ConstAccessor>
maskImage(FileMapImage<PixelType, Alloc> const & img)
{
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser, 
                typename FileMapImage<PixelType, Alloc>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

template <class PixelType, class Alloc>
inline pair< typename FileMapImage<PixelType, Alloc>::const_traverser,
             typename FileMapImage<PixelType, Alloc>::ConstAccessor>
maskImage(FileMapImage<PixelType, Alloc> const & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "maskImage(): ROI rectangle outside image.");
    return pair<typename FileMapImage<PixelType, Alloc>::const_traverser,
                typename FileMapImage<PixelType, Alloc>::ConstAccessor>(img.upperLeft() + ul,
                                                                      img.accessor());
}

} // namespace vigra

#endif