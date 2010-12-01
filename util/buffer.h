#ifndef _BUFFER_H_
#define _BUFFER_H_

#include <vector>
#include <cassert>
#include <utility>

namespace edu_berkeley_cs184 {
namespace util {

template <class T>
class PlaybackBuffer {
public:
	PlaybackBuffer() { _pos = 0; }
	~PlaybackBuffer() {}

	/** Append T to this buffer */
	inline void append(const T&);

	/** gets the current element pointed to by pos */
	inline T getCurrent() const;

	/** sets pos = 0. does NOT reset the buffer */
	inline void resetPointer();

	/** resets both the pointer to 0 and the buffer to be empty */
	inline void resetBuffer();

	/** returns pos() == size() */
	inline bool isEndOfBuffer() const;

	/** advance the pointer. throws an exception if pointer is advanced past the
	 * end of the buffer */
	inline void advance();

	/** returns the current pointer */
	inline size_t pos() const;

	/** returns the size of the buffer */
	inline size_t size() const;

private:
	std::vector<T> _buffer;
	size_t _pos;
};

template <class T>
inline void PlaybackBuffer<T>::append(const T& _t) {
	_buffer.push_back(_t);
}

template <class T>
inline T PlaybackBuffer<T>::getCurrent() const { 
	assert( !isEndOfBuffer() );
	return _buffer[_pos]; 
}

template <class T>
inline void PlaybackBuffer<T>::resetPointer() { _pos = 0; }

template <class T>
inline void PlaybackBuffer<T>::resetBuffer() { _pos = 0; _buffer.clear(); }

template <class T>
inline bool PlaybackBuffer<T>::isEndOfBuffer() const { return _pos == _buffer.size(); }

template <class T>
inline void PlaybackBuffer<T>::advance() { 
	assert( !isEndOfBuffer() );
	_pos++;
}

template <class T>
inline size_t PlaybackBuffer<T>::pos() const { return _pos; }

template <class T>
inline size_t PlaybackBuffer<T>::size() const { return _buffer.size(); }

template <class T0, class T1, class T2>
class Tuple3 {
public:
	Tuple3(const T0& _t0,
				 const T1& _t1,
				 const T2& _t2) : _1(_t0), _2(_t1), _3(_t2) {}
	~Tuple3() {}
	const T0 _1;
	const T1 _2;
	const T2 _3;
};

template <class T0, class T1>
class PlaybackBufferContainer2 {
public:
	PlaybackBufferContainer2() {}
	~PlaybackBufferContainer2() {}

	inline void append(const T0&, const T1&);
	inline std::pair<T0, T1> getCurrent() const;
	inline void resetPointer();
	inline void resetBuffer();
	inline bool isEndOfBuffer() const;
	inline void advance();
	inline size_t pos() const;
	inline size_t size() const;
private:
	PlaybackBuffer<T0> _buf0;
	PlaybackBuffer<T1> _buf1;
};

template <class T0, class T1>
inline void PlaybackBufferContainer2<T0, T1>::append(const T0& _t0, const T1& _t1) {
	_buf0.append(_t0); 
	_buf1.append(_t1);
}

template <class T0, class T1>
inline std::pair<T0, T1> PlaybackBufferContainer2<T0, T1>::getCurrent() const {
	return std::pair<T0, T1>(_buf0.getCurrent(), _buf1.getCurrent());
}

template <class T0, class T1>
inline void PlaybackBufferContainer2<T0, T1>::resetPointer() {
	_buf0.resetPointer();
	_buf1.resetPointer();
}

template <class T0, class T1>
inline void PlaybackBufferContainer2<T0, T1>::resetBuffer() {
	_buf0.resetBuffer();
	_buf1.resetBuffer();
}

template <class T0, class T1>
inline bool PlaybackBufferContainer2<T0, T1>::isEndOfBuffer() const {
	return _buf0.isEndOfBuffer();
}

template <class T0, class T1>
inline void PlaybackBufferContainer2<T0, T1>::advance() {
	_buf0.advance();
	_buf1.advance();
}

template <class T0, class T1>
inline size_t PlaybackBufferContainer2<T0, T1>::pos() const {
	return _buf0.pos();
}

template <class T0, class T1>
inline size_t PlaybackBufferContainer2<T0, T1>::size() const {
	return _buf0.size();
}

template <class T0, class T1, class T2>
class PlaybackBufferContainer3 {
public:
	PlaybackBufferContainer3() {}
	~PlaybackBufferContainer3() {}

	inline void append(const T0&, const T1&, const T2&);
	inline Tuple3<T0, T1, T2> getCurrent() const;
	inline void resetPointer();
	inline void resetBuffer();
	inline bool isEndOfBuffer() const;
	inline void advance();
	inline size_t pos() const;
	inline size_t size() const;
private:
	PlaybackBuffer<T0> _buf0;
	PlaybackBuffer<T1> _buf1;
	PlaybackBuffer<T2> _buf2;
};

template <class T0, class T1, class T2>
inline void PlaybackBufferContainer3<T0, T1, T2>::append(const T0& _t0, const T1& _t1, const T2& _t2) {
	_buf0.append(_t0); 
	_buf1.append(_t1);
	_buf2.append(_t2);
}

template <class T0, class T1, class T2>
inline Tuple3<T0, T1, T2> PlaybackBufferContainer3<T0, T1, T2>::getCurrent() const {
	return Tuple3<T0, T1, T2>(_buf0.getCurrent(), _buf1.getCurrent(), _buf2.getCurrent());
}

template <class T0, class T1, class T2>
inline void PlaybackBufferContainer3<T0, T1, T2>::resetPointer() {
	_buf0.resetPointer();
	_buf1.resetPointer();
	_buf2.resetPointer();
}

template <class T0, class T1, class T2>
inline void PlaybackBufferContainer3<T0, T1, T2>::resetBuffer() {
	_buf0.resetBuffer();
	_buf1.resetBuffer();
	_buf2.resetBuffer();
}

template <class T0, class T1, class T2>
inline bool PlaybackBufferContainer3<T0, T1, T2>::isEndOfBuffer() const {
	return _buf0.isEndOfBuffer();
}

template <class T0, class T1, class T2>
inline void PlaybackBufferContainer3<T0, T1, T2>::advance() {
	_buf0.advance();
	_buf1.advance();
	_buf2.advance();
}

template <class T0, class T1, class T2>
inline size_t PlaybackBufferContainer3<T0, T1, T2>::pos() const {
	return _buf0.pos();
}

template <class T0, class T1, class T2>
inline size_t PlaybackBufferContainer3<T0, T1, T2>::size() const {
	return _buf0.size();
}

}
}

#endif
