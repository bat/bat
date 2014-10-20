#ifndef __IntegrationModel__h
#define __IntegrationModel__h 1

#include <BAT/BCModel.h>

#include <vector>

/// todo documentation, paper reference
class IntegrationModel : public BCModel {
public:
	IntegrationModel();
	~IntegrationModel();

	void SetDimensionality(unsigned int k);
	const unsigned & GetDimensionality() const
	{ return fDimensionality; }

	void SetModality(unsigned int d)
	{ fModality=d; }
	const unsigned & GetModality() const
	{ return fModality; }

	void SetComplexity(unsigned int m)
	{ fComplexity=m; }
	const unsigned & GetComplexity() const
	{ return fComplexity; }

	bool DefineParameters();

	double Phi(unsigned int i, std::vector<double> x) const;
	double S(unsigned int i, double x) const;

	void PopulatePolynomialDegrees();

	double Integral();

	double Likelihood(const std::vector<double> &parameters);
	double LogLikelihood(const std::vector<double> &parameters);

protected:
	unsigned int fDimensionality;
	unsigned int fModality;
	unsigned int fComplexity;

	std::vector<std::vector<unsigned int> > fPolynomialDegrees;
	std::vector<unsigned int> fPolynomialModality;

	std::vector<unsigned int> Divisors(unsigned int d);
	const unsigned & Degree(const unsigned & d) {return d;}
	const unsigned & GetPolynomialDegree(const unsigned & m, const unsigned & k) const
	{ return fPolynomialDegrees[m][k]; }
	std::vector<unsigned int> GetPolynomialDegrees(const unsigned & m) const
   { return fPolynomialDegrees[m]; }
	std::vector<std::vector<unsigned int> > GetPolynomialDegrees() const
	{return fPolynomialDegrees;}

	const unsigned & GetPolynomialModality(const unsigned & m) {return fPolynomialModality[m];}
};

#endif
