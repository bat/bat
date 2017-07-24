#ifndef __BCLOG__H
#define __BCLOG__H

/**
 * @class BCLog
 * @brief A class for managing log messages.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class manages log messages for printing on the screen
 * and into a log file
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

/**
 * Macros to output a message to the log including the class and method name.
 */
#define BCLOG_INTERNAL_OUT(t, s) \
    do {\
    BCLog:: t \
      (std::string(__PRETTY_FUNCTION__) + ": " + s); \
    } while (false)

#define BCLOG_DEBUG(s)   BCLOG_INTERNAL_OUT(OutDebug, s)

#define BCLOG_DETAIL(s)  BCLOG_INTERNAL_OUT(OutDetail, s)

#define BCLOG_ERROR(s)   BCLOG_INTERNAL_OUT(OutError, s)

#define BCLOG_SUMMARY(s) BCLOG_INTERNAL_OUT(OutSummary, s)

#define BCLOG_WARNING(s) BCLOG_INTERNAL_OUT(OutWarning, s)

// ---------------------------------------------------------

#include <fstream>
#include <string>

// ---------------------------------------------------------

class BCLog
{
public:

    // definition of log level

    /**
     * Enumerator for the amount of details to put into the log file */
    enum LogLevel {
        debug,                  ///< Print everything, including debug info
        detail,                 ///< Print all details of operation
        summary,                ///< Print only results summary, warnings, and errors
        warning,                ///< Print only warnings and errors
        error,                  ///< Print only errors
        nothing                 ///< Print nothing
    };

    /** \name Constructor and destructor */
    /** @{ */

    /**
     * Constructor. */
    BCLog();

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
    { fMinimumLogLevelFile = loglevel; };

    /**
     * Sets the minimum log level for screen output.
     * @param loglevel log level */
    static void SetLogLevelScreen(BCLog::LogLevel loglevel)
    { fMinimumLogLevelScreen = loglevel; };

    /**
     * Sets the minimum log level for file and screen output.
     * @param loglevelfile log level for file
     * @param loglevelscreen log level for screen */
    static void SetLogLevel(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen)
    { fMinimumLogLevelFile = loglevelfile; fMinimumLogLevelScreen = loglevelscreen; };

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
    static void OpenLog(const std::string& filename = "log.txt", BCLog::LogLevel loglevelfile = BCLog::debug, BCLog::LogLevel loglevelscreen = BCLog::summary);

    /**
     * @returns true if log file is open or false if not. */
    static bool IsOpen()
    { return fOutputStream.is_open(); }

    /**
     * Closes the log file */
    static void CloseLog()
    { fOutputStream.close(); }

    /**
     * Writes string to the file and screen log if the log level is equal or greater than the minimum
     * @param loglevelfile loglevel for the current message
     * @param loglevelscreen loglevel for the current message
     * @param message string to write to the file and screen log */
    static void Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const std::string& message);

    static void Out(const std::string& message)
    { Out(BCLog::fMinimumLogLevelFile, BCLog::fMinimumLogLevelScreen, message); }

    static void Out(BCLog::LogLevel loglevel, const std::string& message)
    { Out(loglevel, loglevel, message); };

    static void OutError(const std::string& message)
    { Out(error, message); };

    static void OutWarning(const std::string& message)
    { Out(warning, message); };

    static void OutSummary(const std::string& message)
    { Out(summary, message); };

    static void OutDetail(const std::string& message)
    { Out(detail, message); };

    static void OutDebug(const std::string& message)
    { Out(debug, message); };

    /**
     * Writes startup information onto screen and into a logfile */
    static void StartupInfo();

    /**
     * @return string containing the version number  */
    static const std::string& GetVersion()
    { return fVersion; };

    /**
     * Converts a log level to a string */
    static std::string ToString(BCLog::LogLevel);

    /** @} */
private:

    /**
     * BAT version number */
    static std::string fVersion;

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

};

// ---------------------------------------------------------

#endif
