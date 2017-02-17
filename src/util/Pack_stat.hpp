/*
 * Pack_stat.hpp
 *
 *  Created on: Jul 17, 2015
 *      Author: i-bird
 */

#ifndef SRC_PACK_STAT_HPP_
#define SRC_PACK_STAT_HPP_

/*! \brief Unpacking status object
 *
 *
 */
class Unpack_stat
{
	//! offset
	size_t cnt;

public:

	inline Unpack_stat()
	:cnt(0)
	{}

	/*! \brief Increment the offset pointer by off
	 *
	 * \param off
	 *
	 */
	inline void addOffset(size_t off)
	{
		cnt += off;
	}

	/*! \brief Return the actual counter
	 *
	 * \return the counter
	 *
	 */
	inline size_t getOffset()
	{
		return cnt;
	}
};

/*! \brief Packing status object
 *
 *
 */
class Pack_stat
{
	//! marker used to remember some position
	size_t p_mark;

	//! packing offset
	size_t un_ele;

public:


	inline Pack_stat()
	:p_mark(0),un_ele(0)
	{}

	/*! \brief Increment the request pointer
	 *
	 *
	 */
	inline void incReq()
	{
		un_ele++;
	}

	/*! \brief return the actual request for packing
	 *
	 * \return the actual request for packing
	 *
	 */
	inline size_t reqPack()
	{
		return un_ele;
	}

	/*! \brief Mark
	 *
	 *
	 *
	 */
	inline void mark()
	{
		p_mark = un_ele;
	}

	/*! \brief Return the mark
	 *
	 * \return the mark
	 *
	 */
	inline size_t getMark()
	{
		return p_mark;
	}
};

#endif /* SRC_PACK_STAT_HPP_ */
