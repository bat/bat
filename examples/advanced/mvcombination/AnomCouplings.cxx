#include "AnomCouplings.h"

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
AnomCouplings::AnomCouplings() : MVFit()
{
	DefineParameters();
	DefineObservables();
}

// ---------------------------------------------------------
AnomCouplings::~AnomCouplings()
{
}

// ---------------------------------------------------------
void AnomCouplings::DefineParameters()
{
	// add parameters here
	AddParameter("vl", -1.0, 1.0);
	AddParameter("vr", -1.0, 1.0);
	AddParameter("gl", -0.5, 0.5);
	AddParameter("gr", -0.5, 1.0);
	SetPriorDelta("vl", 1.);
	SetPriorDelta("vr", 0.);
	SetPriorConstant("gl");
	SetPriorConstant("gr");
}

// ---------------------------------------------------------
void AnomCouplings::DefineObservables()
{
	AddObservable("F0", 0., 1.);
	AddObservable("FL", 0., 1.);
}

// ---------------------------------------------------------
double AnomCouplings::CalculateObservable(int index, const std::vector<double> &parameters)
{

	double vl = parameters[0];
	double vr = parameters[1];
	double gl = parameters[2];
	double gr = parameters[3];
	
	double xmt = 172.5;
	double xmw = 80.399;
	double xmb = 4.8;
	
	double eb = (xmt*xmt-xmb*xmb-xmw*xmw)/(2.*xmw);
	double qb = sqrt(eb*eb-xmb*xmb);
	
	double xw = xmw/xmt;
	double xb = xmb/xmt;
	
	double q = sqrt(xmt*xmt*xmt*xmt + xmw*xmw*xmw*xmw +  
									xmb*xmb*xmb*xmb - 2.*(xmt*xmw)*(xmt*xmw) -
									2.*(xmt*xmb)*(xmt*xmb) - 2.*(xmb*xmw)*(xmb*xmw))/(2.*xmt);
	
	double Gamma_0 = q*((vl*vl+vr*vr)*(1.-xw*xw-2.*xb*xb-
																		 (xw*xb)*(xw*xb)+xb*xb)/(xw*xw)+
											(gl*gl+gr*gr)*(1.-xw*xw+xb*xb)-
											4.*xb*(vl*vr+gl*gr)-
											2.*vl*(gr-xb*gl)*(1.-xb*xb)/xw-
											2.*vr*(gl-xb*gr)*(1.-xb*xb)/xw+
											2.*xw*vl*(gr+xb*gl)+2.*xw*vr*(gl+xb*gr));
	
	double g1 = q*((vl*vl+vr*vr)*(1.-xw*xw+xb*xb)-
								 4.*xb*(vl*vr+gl*gr)+
								 (gl*gl+gr*gr)*(1.-xw*xw-2.*xb*xb-
																(xw*xb)*(xw*xb)+xb*xb)/(xw*xw)-
								 2.*vl*(gr-xb*gl)*(1.-xb*xb)/xw-
								 2.*vr*(gl-xb*gr)*(1.-xb*xb)/xw+
								 2.*xw*vl*(gr+xb*gl)+2.*xw*vr*(gl+xb*gr));
	
	double g2 = (xmt*xmt*xmt/(2.*xmw*xmw))*(
																					(gl*gl-gr*gr)*(1.-xb*xb)+2.*xw*vl*(gr+xb*gl)-
																					(xw*xw)*(vl*vl-vr*vr)-
																					2.*xw*vr*(gl+xb*gr))*
		(1.-2.*xw*xw-2.*xb*xb-2.*(xb*xw)*(xb*xw)+
		 xw*xw*xw*xw+xb*xb*xb*xb);
	
	double Gamma_p = g1 + g2;
	double Gamma_m = g1 - g2;
	
	double aux = Gamma_0 + Gamma_m + Gamma_p;
	double f0 = Gamma_0/aux;
	double fl = Gamma_m/aux;
	double fr = Gamma_p/aux;
	
	if (index == 0)
		return f0;
	else if (index == 1)
		return fl;
	else return -1;
}

// ---------------------------------------------------------
