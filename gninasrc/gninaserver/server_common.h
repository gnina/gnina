/*
 * server_common.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#ifndef SERVER_COMMON_H_
#define SERVER_COMMON_H_

#include <string>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include <boost/bind.hpp>
#include <boost/assign.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>

#include <fstream>
using namespace std;

typedef boost::shared_ptr<boost::asio::ip::tcp::iostream> stream_ptr;

#endif /* SERVER_COMMON_H_ */
