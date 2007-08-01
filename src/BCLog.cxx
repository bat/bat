#include "BCLog.h" 

std::ofstream BCLog::fOutputStream; 

BCLog::LogLevel BCLog::fMinimumLogLevelFile = BCLog::debug; 

BCLog::LogLevel BCLog::fMinimumLogLevelScreen = BCLog::summary; 

// --------------------------------------------------------- 

BCLog::BCLog()
{

}

// --------------------------------------------------------- 

BCLog::~BCLog()
{

}

// --------------------------------------------------------- 

void BCLog::OpenLog(const char* filename, BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen)
{

  // open log file 

  BCLog::fOutputStream.open(filename); 

  if (!BCLog::fOutputStream.is_open())
    {
      std::cout << " Could not open log file " << filename << ". " << std::endl; 
      return; 
    }

  // set log level 

  BCLog::SetMinimumLogLevelFile(loglevelfile); 

  BCLog::SetMinimumLogLevelScreen(loglevelscreen); 

}

// --------------------------------------------------------- 

void BCLog::OpenLog(const char* filename)
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

void BCLog::Out(BCLog::LogLevel loglevelfile, BCLog::LogLevel loglevelscreen, const char* message)
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

void BCLog::Out(const char* message)
{

  BCLog::Out(BCLog::fMinimumLogLevelFile, BCLog::fMinimumLogLevelScreen, message); 

}

// --------------------------------------------------------- 

char* BCLog::ToString(BCLog::LogLevel loglevel)
{

  if (loglevel == debug) 
    return "Debug"; 

  else if (loglevel == 1) 
    return "Detail"; 

  else if (loglevel == 2) 
    return "Summary"; 

  else if (loglevel == 3) 
    return "Warning"; 

  else
    return ""; 

}

// --------------------------------------------------------- 

