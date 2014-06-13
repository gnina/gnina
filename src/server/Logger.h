/*
 * Logger.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 *
 *  A multi-threading aware class for logging concurrently to a file.
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

class Logger
{
	FILE *LOG;
	boost::mutex locker;
public:
	Logger(const string& logname): LOG(NULL)
	{
		if(logname.size() == 0) return;
		LOG = fopen(logname.c_str(), "a");
		if (LOG == NULL)
		{
			cerr << "Could not open log directory " << logname << "\n";
			exit(-1);
		}
		setlinebuf(LOG);
	}

	~Logger()
	{
		fclose(LOG);
	}

	//output message with printfstyle arguments atomically to log
	void log(const char* str, ...)
	{
		if(!LOG) return;
		va_list argptr;
		va_start(argptr, str);

		locker.lock();

		vfprintf(LOG, str, argptr);

		locker.unlock();
		va_end(argptr);
	}
};

#endif /* LOGGER_H_ */
