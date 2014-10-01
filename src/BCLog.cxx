/*
 * Copyright (C) 2007-2014, the BAT core developer team
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

#include <fstream>
#include <iomanip>


std::ofstream BCLog::fOutputStream;

BCLog::LogLevel BCLog::fMinimumLogLevelFile = BCLog::debug;

BCLog::LogLevel BCLog::fMinimumLogLevelScreen = BCLog::summary;

bool BCLog::fFirstOutputDone = false;

const char * BCLog::fVersion = VERSION;

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
      std::cerr << " Could not open log file " << filename << ". " << std::endl;
      return;
   }

   // set log level
   BCLog::SetLogLevelFile(loglevelfile);
   BCLog::SetLogLevelScreen(loglevelscreen);

   BCLog::Out(BCLog::summary,BCLog::summary,Form("Opening logfile %s",filename));
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

bool BCLog::IsOpen()
{
   return BCLog::fOutputStream.is_open();
}

// ---------------------------------------------------------

void BCLog::CloseLog()
{
   BCLog::fOutputStream.close();
}

// ---------------------------------------------------------

void BCLog::Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const char * message)
{
   // if this is the first call to Out(), call StartupInfo() first
   if(!fFirstOutputDone)
      BCLog::StartupInfo();

   // open log file if not opened
   if (BCLog::IsOpen())
   {
      // write message in to log file
      if (loglevelfile >= BCLog::fMinimumLogLevelFile)
         BCLog::fOutputStream << BCLog::ToString(loglevelfile) << " : " << message << std::endl;
   }

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
   char * message = Form(
                         " +------------------------------------------------------+\n"
                         " |                                                      |\n"
                         " | BAT version %7s                                  |\n"
                         " | Copyright (C) 2007-2014, the BAT core developer team |\n"
                         " | All rights reserved.                                 |\n"
                         " |                                                      |\n"
                         " | For the licensing terms see doc/COPYING              |\n"
                         " | For documentation see http://mpp.mpg.de/bat          |\n"
                         " |                                                      |\n"
                         " +------------------------------------------------------+\n",
                         BCLog::fVersion);
   
   // write message to screen
   if (BCLog::fMinimumLogLevelScreen < BCLog::nothing)
      std::cout << message << std::endl;

   if (BCLog::IsOpen() && BCLog::fMinimumLogLevelFile<BCLog::nothing)
      BCLog::fOutputStream << message;

   fFirstOutputDone = true;
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
      case error:
         return "Error  ";
      default:
         return "";
   }
}

// ---------------------------------------------------------

int printBATUponLoading()
{
  /*
    std::cout <<
    " +------------------------------------------------------+\n"
    " |                                                      |\n"
    " | BAT version " << std::setw(7) << VERSION << "                                  |\n"
    " | Copyright (C) 2007-2014, the BAT core developer team |\n"
    " | All rights reserved.                                 |\n"
    " |                                                      |\n"
    " | For the licensing terms see doc/COPYING              |\n"
    " | For documentation see http://mpp.mpg.de/bat          |\n"
    " |                                                      |\n"
    " +------------------------------------------------------+\n";
  */
   return 0;
}

static int tmpvarPrint = printBATUponLoading();
