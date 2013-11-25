#pragma once
#include <vector>

//-----------------------------------------------------------------------------
//! The dump stream allows a class to record its internal state to a memory object
//! so that it can be restored later.
//! In FEBio this is used for storing the FEModel state during running restarts
class DumpStream
{
public:
	DumpStream();
	~DumpStream();

	void clear();

	void set_position(int l);

	void write(void* pd, int nsize);
	void read (void* pd, int nsize);

	template <typename T> DumpStream& operator << (T& o);
	template <typename T> DumpStream& operator >> (T& o);

	template <typename T> DumpStream& operator << (std::vector<T>& o);
	template <typename T> DumpStream& operator >> (std::vector<T>& o);

protected:
	void grow_buffer(int l);

private:
	char*	m_pb;			//!< pointer to buffer
	char*	m_pd;			//!< position to insert a new value
	int		m_nsize;		//!< size of stream
	int		m_nreserved;	//!< size of reserved buffer
};

template <typename T> inline DumpStream& DumpStream::operator << (T& o)
{
	int l = sizeof(T);
	write(&o, l);
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator >> (T& o)
{
	int l = sizeof(T);
	read(&o, l);
	return *this;
}

template <typename T> DumpStream& DumpStream::operator << (std::vector<T>& o)
{
	DumpStream& This = *this;
	int N = (int) o.size();
	for (int i=0; i<N; ++i) This << o[i];
	return This;
}

template <typename T> DumpStream& DumpStream::operator >> (std::vector<T>& o)
{
	DumpStream& This = *this;
	int N = (int) o.size();
	for (int i=0; i<N; ++i) This >> o[i];
	return This;
}
