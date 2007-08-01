#include "BCDataPoint.h" 

#include "TROOT.h"

// --------------------------------------------------------- 

BCDataPoint::BCDataPoint(int nvariables)
{

  fData.assign(nvariables, 0.0); 

}

// --------------------------------------------------------- 

BCDataPoint::BCDataPoint(vector<double> x)
{

  int nentries = int(x.size()); 

  for (int i = 0; i < nentries; i++)
    {
      fData.push_back(x.at(i)); 
    }

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

}

// --------------------------------------------------------- 

void BCDataPoint::SetValues(std::vector <double> values) 
{

  if (values.size() == fData.size())
    for (int i = 0; i < int(fData.size()); i++)
      fData[i] = values[i]; 

}

// --------------------------------------------------------- 


