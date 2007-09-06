/*! \class BCLog
 *  \brief Log printouts
 *
 * This class manages printouts on the screen and into a log file 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 16.05.2007 by Kevin 
 * 
 * REVISION: 
 *
 * --------------------------------------------------------- 
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCLOG__H
#define __BCLOG__H

#include <iostream>
#include <fstream> 

#include <TROOT.h>
#include <TError.h>

#include "BCReleaseVersion.h"

// --------------------------------------------------------- 

class BCLog
{

 public:

	// definition of log level 

	/**
	 * Enumerator for the amount of details to put into the log file 
	 * Log levels: 
	 * debug   : Lowest level of information 
	 * detail  : Details of functions, etc. 
	 * summary : Results 
	 * warning : Warning messages 
	 * nothing : Not written into file 
	 */ 
	enum LogLevel {debug, detail, summary, warning, nothing}; 

	// constructor and destructor 

	/**
	 * A constructor. 
	 */ 
	BCLog(); 

	/**
	 * A destructor. 
	 */ 
	~BCLog(); 

	// methods (get) 

	/** 
	 * Returns the minimum log level for file output. 
	 * @return The log level 
	 */ 
	static BCLog::LogLevel GetMinimumLogLevelFile()
	{ return fMinimumLogLevelFile; }; 

	/** 
	 * Returns the minimum log level for screen output. 
	 * @return The log level 
	 */ 
	static BCLog::LogLevel GetMinimumLogLevelScreen()
	{ return fMinimumLogLevelScreen; }; 

	// method (set) 

	/** 
	 * Sets the minimum log level for file output. 
	 * @param loglevel The log level 
	 */ 
	static void SetMinimumLogLevelFile(BCLog::LogLevel loglevel)
	{ BCLog::fMinimumLogLevelFile = loglevel; }; 

	/** 
	 * Sets the minimum log level for screen output. 
	 * @param loglevel The log level 
	 */ 
	static void SetMinimumLogLevelScreen(BCLog::LogLevel loglevel)
	{ BCLog::fMinimumLogLevelScreen = loglevel; }; 

	// methods 

	/**
	 * Opens log file and sets minimum log levels for file and screen output. 
	 * @param filename The log filename 
	 * @param loglevelfile The minimum log level for file output
	 * @param loglevelscreen The minimum log level for screen output
	 */ 
	static void OpenLog(const char* filename, BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen); 

	static void OpenLog(const char* filename); 

	static void OpenLog(); 

	/**
	 * Closes the log file 
	 */ 
	static void CloseLog(); 

	/*
	 * Writes string to the file and screen log if the log level is equal or greater than the minimum
	 * @param loglevelfile The loglevel for the current message 
	 * @param loglevelscreen The loglevel for the current message 
	 * @param message The string to write to the file and screen log 
	 */ 
	static void Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const char* message); 

	static void Out(const char* message); 

	/**
	 * Writes startup information onto screen and into a logfile
	 */ 
	static void StartupInfo();

	/**
	 * @return string containing the version number
	 */
	static char * GetVersion()
	{ return fVersion; };

	/**
	 * @return unique number for use in histogram name string
	 */
	static int GetHIndex()
	{ return BCLog::fHindex++; };

 private: 

	/**
	 * BAT version number
	 */
	static char * fVersion;

	/**
	 * The minimum file log level 
	 */ 
	static BCLog::LogLevel fMinimumLogLevelFile; 

	/** 
	 * The minimum screen log level
	 */ 
	static BCLog::LogLevel fMinimumLogLevelScreen; 

	/** 
	 * The output stream for the file log 
	 */ 
	static std::ofstream fOutputStream; 

	/**
	 * Converts a log level to a string
	 */ 
	static char* ToString(BCLog::LogLevel); 

	/**
	 * Global histogram counter
	 */ 
	static int fHindex;

}; 

// --------------------------------------------------------- 

#endif 
