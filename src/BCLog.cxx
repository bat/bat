/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "config.h"

#include "BCLog.h"

#include <TError.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>


std::ofstream BCLog::fOutputStream;

BCLog::LogLevel BCLog::fMinimumLogLevelFile = BCLog::debug;

BCLog::LogLevel BCLog::fMinimumLogLevelScreen = BCLog::summary;

bool BCLog::fFirstOutputDone = false;

bool BCLog::fPrefix = true;

std::string BCLog::fVersion = VERSION;

// ---------------------------------------------------------

BCLog::BCLog()
{
    // suppress the ROOT Info printouts
    gErrorIgnoreLevel = 2000;
}


// ---------------------------------------------------------

void BCLog::OpenLog(const std::string& filename, BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen)
{
    // suppress the ROOT Info printouts
    gErrorIgnoreLevel = 2000;

    // first close and flush and existing log file
    BCLog::CloseLog();

    // open log file
    BCLog::fOutputStream.open(filename.data());

    if (!BCLog::fOutputStream.is_open()) {
        std::cerr << " Could not open log file " << filename << ". " << std::endl;
        return;
    }

    // set log level
    BCLog::SetLogLevelFile(loglevelfile);
    BCLog::SetLogLevelScreen(loglevelscreen);

    BCLog::Out(BCLog::summary, BCLog::summary, "Opening logfile " + filename);
}

// ---------------------------------------------------------

void BCLog::Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const std::string& message)
{
    // if this is the first call to Out(), call StartupInfo() first
    if (!fFirstOutputDone)
        BCLog::StartupInfo();

    // open log file if not opened
    if (BCLog::IsOpen()) {
        // write message in to log file
        if (loglevelfile >= BCLog::fMinimumLogLevelFile)
            BCLog::fOutputStream << BCLog::ToString(loglevelfile) << message << std::endl;
    }

    // write message to screen
    if (loglevelscreen >= BCLog::fMinimumLogLevelScreen)
        std::cout << BCLog::ToString(loglevelscreen) << message << std::endl;
}

// ---------------------------------------------------------

void BCLog::StartupInfo()
{
    const char* message = Form(
                              " +------------------------------------------------------+\n"
                              " |                                                      |\n"
                              " | BAT version %-12s                             |\n"
                              " | Copyright (C) 2007-2018, the BAT core developer team |\n"
                              " | All rights reserved.                                 |\n"
                              " |                                                      |\n"
                              " | For the licensing terms see doc/COPYING              |\n"
                              " | For documentation see http://mpp.mpg.de/bat          |\n"
                              " | Please cite: DOI 10.1016/j.cpc.2009.06.026           |\n"
                              " |              http://arxiv.org/abs/0808.2552          |\n"
                              " |                                                      |\n"
                              " +------------------------------------------------------+\n",
                              BCLog::fVersion.data());

    // write message to screen
    if (BCLog::fMinimumLogLevelScreen < BCLog::nothing)
        std::cout << message << std::endl;

    if (BCLog::IsOpen() && BCLog::fMinimumLogLevelFile < BCLog::nothing)
        BCLog::fOutputStream << message;

    fFirstOutputDone = true;
}

// ---------------------------------------------------------

std::string BCLog::ToString(BCLog::LogLevel loglevel)
{
    if (!fPrefix)
        return "";

    switch (loglevel) {
        case debug:
            return "Debug   : ";
        case detail:
            return "Detail  : ";
        case summary:
            return "Summary : ";
        case warning:
            return "Warning : ";
        case error:
            return "Error   : ";
        default:
            return "";
    }
}
