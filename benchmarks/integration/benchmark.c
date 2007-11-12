#include <BCModelGauss.h>
#include <BCModelManager.h> 
#include <BCLog.h>

#include <BCTest.h> 

#include <TStopwatch.h> 

// ---------------------------------------------------------
  
int main()
{

/* 	BCParameter * para1 = new BCParameter();  */

/* 	cout << " Parameter 1 : " << endl;  */
/* 	cout << " Name        : " << para1 -> GetName() << endl;  */
/* 	cout << " Lower limit : " << para1 -> GetLowerLimit() << endl;  */
/* 	cout << " Upper limit : " << para1 -> GetUpperLimit() << endl;  */
	
/* 	delete para1;  */

	BCModel * model1 = new BCModel(); 
	BCModel * model2 = new BCModel(); 

	BCModelManager * manager1 = new BCModelManager; 
	manager1 -> AddModel(model1); 
	manager1 -> AddModel(model2); 

	BCModelManager * manager2(manager1); 

	cout << manager1 -> GetNModels() << endl; 
	cout << manager2 -> GetNModels() << endl; 
	
//	BCTest * testy = new BCTest("test"); 
//	delete testy; 

//	BCModel * model = new BCModel(); 
//	delete model; 

	return 0; 
}

/* 	// --------------------------------------------------------- */
/* 	// create stopwatch  */
/* 	// --------------------------------------------------------- */

/* 	TStopwatch * fStopwatch = new TStopwatch();  */

/* 	// --------------------------------------------------------- */
/* 	// open log file  */
/* 	// --------------------------------------------------------- */

/* 	BCLog::OpenLog();  */

/* 	// --------------------------------------------------------- */
/* 	// model definition  */
/* 	// --------------------------------------------------------- */

/* 	BCModelGauss* fModelGauss = new BCModelGauss("ModelGauss");  */

/* 	BCParameter * para = fModelGauss -> GetParameter(0);  */

/* 	BCParameter * para1 = new BCParameter();  */
	
/* 	delete para1;  */

/* 	return 0;  */
/* 	//	BCParameter * para1 =fModelGauss -> GetParameter(0);  */

/* 	if (para == para1)  */
/* 		cout << "same" <<endl;  */


/* 	return 1;  */

/* 	// --------------------------------------------------------- */
/* 	// do normalization with different methods  */
/* 	// --------------------------------------------------------- */

/* 	// Monte Carlo integration  */

/* 	// start stopwatch  */

/* 	std::cout << std::endl;  */
/* 	std::cout << " Integrate with Monte Carlo method " << endl;  */

/* 	fModelGauss -> SetIntegrationMethod(BCIntegrate::kIMonteCarlo);  */
/* 	fModelGauss -> SetNIterationsMax(10000000);  */
/* 	fModelGauss -> SetRelativePrecision(1e-2);  */

/* 	fStopwatch -> Start();  */

/* 	double normalization_montecarlo = fModelGauss -> Normalize();  */

/* 	fStopwatch -> Stop();  */

/* 	double realtime_montecarlo = fStopwatch -> RealTime();  */
/* 	double cputime_montecarlo  = fStopwatch -> CpuTime();  */

/* 	// Monte Carlo integration  */

/* 	// start stopwatch  */

/* 	std::cout << std::endl;  */
/* 	std::cout << " Integrate with method from Cuba library (VEGAS) " << endl;  */

/* 	fModelGauss -> SetIntegrationMethod(BCIntegrate::kICuba);  */
/* 	fModelGauss -> SetNIterationsMax(10000000);  */
/* 	fModelGauss -> SetRelativePrecision(1e-2);  */

/* 	fStopwatch -> Start(kTRUE);  */

/* 	double normalization_vegas = fModelGauss -> Normalize();  */

/* 	cout << " norm " << normalization_vegas << endl;  */

/* 	fStopwatch -> Stop();  */

/* 	double realtime_vegas = fStopwatch -> RealTime();  */
/* 	double cputime_vegas  = fStopwatch -> CpuTime();  */

/* 	// other integration methods ...  */

/* 	// ..  */

/* 	// --------------------------------------------------------- */
/* 	// print results of benchmark test */
/* 	// --------------------------------------------------------- */

/* 	std::cout << std::endl;  */
/* 	std::cout << " Summary of integration methods " << std::endl;  */
/* 	std::cout << " ============================== " << std::endl;  */
/* 	std::cout << std::endl;  */

/* 	std::cout << " Monte Carlo " << std::endl;  */
/* 	std::cout << " Integral   :  " << normalization_montecarlo << std::endl; */
/* 	std::cout << " Real time  :  " << realtime_montecarlo << std::endl;  */
/* 	std::cout << " CPU time   :  " << cputime_montecarlo << std::endl;  */
/* 	std::cout << std::endl;  */

/* 	std::cout << " VEGAS " << std::endl;  */
/* 	std::cout << " Integral   :  " << normalization_vegas << std::endl; */
/* 	std::cout << " Real time  :  " << realtime_vegas << std::endl;  */
/* 	std::cout << " CPU time   :  " << cputime_vegas << std::endl;  */
/* 	std::cout << std::endl;  */

/* 	// --------------------------------------------------------- */
/* 	// close log file  */
/* 	// --------------------------------------------------------- */

/* 	// closes the log file  */

/* 	BCLog::CloseLog();  */

/* 	return 0;  */

/* } */

// ---------------------------------------------------------
  
