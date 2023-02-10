/*
 * ObjReader.hpp
 *
 *  Created on: Oct 31, 2018
 *      Author: i-bird
 */

#ifndef OBJREADER_HPP_
#define OBJREADER_HPP_

#include "config.h"
#include "Space/Shape/Point.hpp"

#ifdef HAVE_TINYOBJLOADER
#include "tiny_obj_loader.h"

/*! \brief Wavefront obj File reader
 *
 * \tparam T precision of the 3D space
 *
 */
template<typename T>
class ObjReader
{
	size_t shape_counter;
	size_t face_counter;
	size_t index_offset;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	size_t end_shape;

	public:

	/*! \brief Constructor
	 *
	 *
	 */
	ObjReader()
	{}

	/*! \brief Read an Wavefront obj file
	 *
	 * \param file to open
	 *
	 * \return true if the read succeed
	 */
	bool read(std::string file)
	{
		shape_counter = 0;
		face_counter = 0;
		index_offset = 0;

		shapes.clear();
		materials.clear();

		std::string err;
		bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, file.c_str());

		if (!err.empty())
		{std::cerr << err << std::endl;}

		if (!ret)
		{return false;}

		end_shape = shapes.size();

		return true;
	}


	/*! \brief Next triangle
	 *
	 * \return itself
	 *
	 */
	ObjReader & operator++()
	{
		index_offset += getFaceNVertex();
		++face_counter;

		if (face_counter >= shapes[shape_counter].mesh.num_face_vertices.size())
		{
			face_counter = 0;
			index_offset = 0;
			++shape_counter;
		}

		return *this;
	}

	/*! \brief Return true if we have another face
	 *
	 * \return true if we have next face
	 */
	inline bool isNext()
	{
		return shape_counter < end_shape;
	}

	/*! \brief Iterate only on a partucular object faces
	 *
	 * \param i object to iterate
	 *
	 * \return
	 */
	inline void setObject(int i)
	{
		if (i >= shapes.size())
		{
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " you selected object " << i << " but the file contain " << shapes.size() << " objects" << std::endl;
		}

		shape_counter = i;
		end_shape = i+1;
	}


	/*! \brief Get the number of vertices in the face
	 *
	 * \return number of vertices a face has
	 */
	unsigned int getFaceNVertex()
	{
		return  shapes[shape_counter].mesh.num_face_vertices[face_counter];
	}

	/*! \brief return the vertex v of the actual face
	 *
	 * \param v vertex
	 * \return the position of such vertex
	 */
	inline Point<3,T> getVertex(unsigned int v)
	{
		Point<3,T> p;

		tinyobj::index_t idx = shapes[shape_counter].mesh.indices[index_offset + v];

		p.get(0) = attrib.vertices[3*idx.vertex_index+0];
		p.get(1) = attrib.vertices[3*idx.vertex_index+1];
		p.get(2) = attrib.vertices[3*idx.vertex_index+2];

		return p;
	}

	/*! \brief return the vertex v of the actual face
	 *
	 * \param v vertex
	 * \return the position of such vertex
	 */
	inline Point<3,T> getNormal(unsigned int v)
	{
		Point<3,T> p;

		tinyobj::index_t idx = shapes[shape_counter].mesh.indices[index_offset + v];

		p.get(0) = attrib.normals[3*idx.normal_index+0];
		p.get(1) = attrib.normals[3*idx.normal_index+1];
		p.get(2) = attrib.normals[3*idx.normal_index+2];

		return p;
	}

	/*! \brief Get the texture coordinates
	 *
	 * \param v vertex
	 *
	 * \return The texture coordinates of the vertex v
	 *
	 */
	inline Point<2,T> getTexCoord(unsigned int v)
	{
		Point<2,T> p;

		tinyobj::index_t idx = shapes[shape_counter].mesh.indices[index_offset + v];

		p.get(0) = attrib.texcoords[2*idx.texcoord_index+0];
		p.get(1) = attrib.texcoords[2*idx.texcoord_index+1];

		return p;
	}
};

#endif

#endif /* OBJREADER_HPP_ */
