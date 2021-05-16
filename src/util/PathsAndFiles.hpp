//
// Created by jstark on 2019-12-03.
//
/**
 * @file PathsAndFiles.hpp
 *
 * @brief Header file containing functions for creating files and folders.
 *
 * @author Justina Stark
 * @date December 2019
 */

#ifndef FILES_READ_AND_WRITE_PATHSANDFILES_HPP
#define FILES_READ_AND_WRITE_PATHSANDFILES_HPP
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "VCluster/VCluster.hpp"
/**@brief Gets the current working directory and returns path as string.
 *
 * @return Std::string of path to current working directory.
 */
std::string get_cwd()
{
	char *cwd = nullptr;
	size_t size;
	cwd = getcwd(cwd, size);
	
	// convert current directory to string
	std::string s_cwd;
	s_cwd = cwd;

	//	std::cout << "The current working directory is: " << dir << std::endl;
	return s_cwd;
}
/**@brief Checks if a file already exists.
 *
 * @param path Std::string with path of file for which existence should be checked.
 * @return True, if file exists, false if not.
 */
bool check_if_file_exists(std::string path)
{
	try
	{
		// Create a filesystem::path object from given path string
		boost::filesystem::path pathObj(path);
		if (boost::filesystem::exists(pathObj) && boost::filesystem::is_regular_file(pathObj))
		{
			return true;
//			BOOST_LOG_TRIVIAL(info) << path << " -> File exists.";
		}
		
	}
	catch (boost::filesystem::filesystem_error & e)
	{
		std::cerr << e.what() << std::endl;
//		BOOST_LOG_TRIVIAL(error) << "Error when checking existence of file ( " << path << " ): " << e.what();
	}
	return false;
}
/**@brief Creates a file if not already existent.
 *
 * @param path Std::string that contains path including filename of the file that should be created.
 */
void create_file_if_not_exist(std::string path)
{
	auto & v_cl = create_vcluster();
	if (v_cl.rank() == 0)
	{
		if ( ! check_if_file_exists(path))
		{
			std::ofstream outfile (path);
//			BOOST_LOG_TRIVIAL(info) << "Created file with name: " << path;
		}
	}
}


/**@brief Checks if a directory already exists.
 *
 * @param path Std::string with path of directory for which existence should be checked.
 * @return True, if directory exists, false if not.
 */
bool check_if_directory_exists(std::string path)
{
	try
	{
		// Create a filesystem::path object from given path string
		boost::filesystem::path pathObj(path);
		if (boost::filesystem::exists(pathObj) && boost::filesystem::is_directory(pathObj)) return true;
	}
	catch (boost::filesystem::filesystem_error & e)
	{
		std::cerr << e.what() << std::endl;
	}
	return false;
}
/**@brief Creates a directory if not already existent.
 *
 * @param path Std::string that contains path including name of the folder that should be created.
 */
void create_directory_if_not_exist(std::string path)
{
	auto & v_cl = create_vcluster();
	if (v_cl.rank() == 0)
	{
		if ( ! check_if_directory_exists(path))
		{
			boost::filesystem::create_directory(path);
//			BOOST_LOG_TRIVIAL(info) << "Created directory with name: " << path;
		}
		else std::cout << "Folder for current settings ( '" << path << "' ) already exists. New files will be saved to this folder." << std::endl;
//		BOOST_LOG_TRIVIAL(info) << "Folder for current settings ( '" << path << "' ) already exists. New files will be saved to this folder.";
	}
}





#endif //FILES_READ_AND_WRITE_PATHSANDFILES_HPP

