#include "BCDataPoint.h" 
#include "BCLog.h" 

#include "TROOT.h"

// --------------------------------------------------------- 

BCDataPoint::BCDataPoint(int nvariables)
{

  fData.assign(nvariables, 0.0); 

}

// --------------------------------------------------------- 

BCDataPoint::BCDataPoint(vector<double> x)
{

  for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); ++it)
    fData.push_back(*it); 

}

// --------------------------------------------------------- 

BCDataPoint::~BCDataPoint()
{

}

// --------------------------------------------------------- 

double BCDataPoint::GetValue(int index) 
{

  double value = -1.0; 

  if (index >= 0 && index < int(fData.size()))
    value = fData[index]; 

  return value; 

}

// --------------------------------------------------------- 

void BCDataPoint::SetValue(int index, double value) 
{

  if (index >= 0 && index < int(fData.size()))
    fData[index] = value; 

  else
    BCLog::Out(BCLog::warning, BCLog::warning,"BCDataPoint::SetValue. Index out of range.");

}

// --------------------------------------------------------- 

void BCDataPoint::SetValues(std::vector <double> values) 
{

  if (values.size() == fData.size())
    for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it)
      fData.push_back(*it); 
  
  else
    BCLog::Out(BCLog::warning, BCLog::warning,"BCDataPoint::SetValues. vectors have different ranges.");
  
}

// --------------------------------------------------------- 


