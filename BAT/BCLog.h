#ifndef __BCLOG__H
#define __BCLOG__H

/*!
 * \class BCLog
 * \brief A class for managing log messages.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class manages log messages for printing on the screen
 * and into a log file
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <iostream>

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
       * error   : Error message
       * nothing : No output
       */
      enum LogLevel {debug, detail, summary, warning, error, nothing};

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * Constructor. */
      BCLog();

      /**
       * Destructor. */
      ~BCLog();

      /** @} */
      /** \name Getters */
      /** @{ */

      /**
       * Returns the minimum log level for file output.
       * @return log level */
      static BCLog::LogLevel GetLogLevelFile()
         { return fMinimumLogLevelFile; };

      /**
       * Returns the minimum log level for screen output.
       * @return log level */
      static BCLog::LogLevel GetLogLevelScreen()
         { return fMinimumLogLevelScreen; };

      /** @} */
      /** \name Setters */
      /** @{ */

      /**
       * Sets the minimum log level for file output.
       * @param loglevel log level */
      static void SetLogLevelFile(BCLog::LogLevel loglevel)
         { BCLog::fMinimumLogLevelFile = loglevel; };

      /**
       * Sets the minimum log level for screen output.
       * @param loglevel log level */
      static void SetLogLevelScreen(BCLog::LogLevel loglevel)
         { BCLog::fMinimumLogLevelScreen = loglevel; };

      /**
       * Sets the minimum log level for file and screen output.
       * @param loglevelscreen log level for screen
       * @param loglevelfile log level for file */
      static void SetLogLevel(BCLog::LogLevel loglevelscreen, BCLog::LogLevel loglevelfile)
         { BCLog::fMinimumLogLevelFile = loglevelfile; BCLog::fMinimumLogLevelScreen = loglevelscreen; };

      /**
       * Sets the minimum log level for file and screen output.
       * @param loglevel log level */
      static void SetLogLevel(BCLog::LogLevel loglevel)
         { SetLogLevel(loglevel, loglevel); };

      /** @} */
      /** \name Miscellaneous */
      /** @{ */

      /**
       * Opens log file and sets minimum log levels for file and screen output.
       * @param filename log filename
       * @param loglevelfile minimum log level for file output
       * @param loglevelscreen minimum log level for screen output */
      static void OpenLog(const char * filename, BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen);

      static void OpenLog(const char * filename);

      static void OpenLog();

      /**
       * @returns true if log file is open or false if not. */
      static bool IsOpen();

      /**
       * Closes the log file */
      static void CloseLog();

      /**
       * Writes string to the file and screen log if the log level is equal or greater than the minimum
       * @param loglevelfile loglevel for the current message
       * @param loglevelscreen loglevel for the current message
       * @param message string to write to the file and screen log */
      static void Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const char * message);

      static void Out(const char * message);

      static void Out(BCLog::LogLevel loglevel, const char * message)
         { Out(loglevel,loglevel,message); };

      static void OutError(const char * message)
         { Out(error,message); };

      static void OutWarning(const char * message)
         { Out(warning,message); };

      static void OutSummary(const char * message)
         { Out(summary,message); };

      static void OutDetail(const char * message)
         { Out(detail,message); };

      static void OutDebug(const char * message)
         { Out(debug,message); };

      /**
       * Writes startup information onto screen and into a logfile */
      static void StartupInfo();

      /**
       * @return string containing the version number  */
      static const char * GetVersion()
         { return fVersion; };

      /**
       * @return unique number for use in histogram name string */
      static int GetHIndex()
         { return BCLog::fHindex++; };

      /**
       * Converts a log level to a string */
      static const char * ToString(BCLog::LogLevel);

      /** @} */
   private:

      /**
       * BAT version number */
      static const char * fVersion;

      /**
       * The minimum file log level */
      static BCLog::LogLevel fMinimumLogLevelFile;

      /**
       * The minimum screen log level */
      static BCLog::LogLevel fMinimumLogLevelScreen;

      /**
       * The output stream for the file log */
      static std::ofstream fOutputStream;

      /**
       * Specifies wheather there were output printouts already */
      static bool fFirstOutputDone;

      /**
       * Global histogram counter */
      static int fHindex;
};

// ---------------------------------------------------------

#endif
