/*   
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.   
 * All rights reserved.   
 *   
 * For the licensing terms see doc/LICENSE.   
 */   
 
// ---------------------------------------------------------  

#include "BCLog.h"

std::ofstream BCLog::fOutputStream;

BCLog::LogLevel BCLog::fMinimumLogLevelFile = BCLog::debug;

BCLog::LogLevel BCLog::fMinimumLogLevelScreen = BCLog::summary;

const char * BCLog::fVersion = BAT_VERSION;

int BCLog::fHindex = 0;

// ---------------------------------------------------------

BCLog::BCLog()
{

	// suppress the ROOT Info printouts

	gErrorIgnoreLevel=2000;

}

// ---------------------------------------------------------

BCLog::~BCLog()
{}

// ---------------------------------------------------------

void BCLog::OpenLog(const char * filename, BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen)
{

	// suppress the ROOT Info printouts

	gErrorIgnoreLevel=2000;

	// open log file

	BCLog::fOutputStream.open(filename);

	if (!BCLog::fOutputStream.is_open())
		{
			std::cout << " Could not open log file " << filename << ". " << std::endl;
			return;
		}

	BCLog::StartupInfo();

	BCLog::Out(BCLog::summary,BCLog::summary,Form("Opening logfile %s",filename));

	// set log level 

	BCLog::SetMinimumLogLevelFile(loglevelfile); 
	BCLog::SetMinimumLogLevelScreen(loglevelscreen); 

}

// --------------------------------------------------------- 

void BCLog::OpenLog(const char * filename)
{

	BCLog::OpenLog(filename, BCLog::debug, BCLog::summary); 

}

// --------------------------------------------------------- 

void BCLog::OpenLog()
{

	BCLog::OpenLog("log.txt", BCLog::debug, BCLog::summary); 

}

// --------------------------------------------------------- 

void BCLog::CloseLog()
{

	BCLog::fOutputStream.close(); 

}

// --------------------------------------------------------- 

void BCLog::Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const char * message)
{

	// open log file if not opened

	if (!BCLog::fOutputStream.is_open())
		{
			std::cout << " Log file not opended. " << std::endl;
			return;
		}

	// write message in to log file

	if (loglevelfile >= BCLog::fMinimumLogLevelFile)
		BCLog::fOutputStream << BCLog::ToString(loglevelfile) << " : " << message << std::endl;

	// write message to screen

	if (loglevelscreen >= BCLog::fMinimumLogLevelScreen)
		std::cout << BCLog::ToString(loglevelscreen) << " : " << message << std::endl;

}

// --------------------------------------------------------- 

void BCLog::Out(const char * message)
{

	BCLog::Out(BCLog::fMinimumLogLevelFile, BCLog::fMinimumLogLevelScreen, message);

}

// --------------------------------------------------------- 

void BCLog::StartupInfo()
{

	// open log file if not opened 

	if (!BCLog::fOutputStream.is_open())
		{
			std::cout << " Log file not opended. " << std::endl; 
			return; 
		}

	char * message = Form(
												" +------------------------------\n"
												" |\n"
												" |     Running with BAT\n"
												" |      Version %s\n"
												" |\n"
												" | http://www.mppmu.mpg.de/bat\n"
												" +------------------------------\n",
												BCLog::fVersion);

	BCLog::fOutputStream << message;
	std::cout << message;

}

// --------------------------------------------------------- 

const char * BCLog::ToString(BCLog::LogLevel loglevel)
{
	
	switch (loglevel)
		{
		case debug:
			return "Debug  ";
		case detail:
			return "Detail ";
		case summary:
			return "Summary";
		case warning:
			return "Warning";
		default:
			return "";
	}

}

// --------------------------------------------------------- 

