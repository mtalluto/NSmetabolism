// // class NSModel {
// // public:
// // 	int getTime(int t);
// // 	void run();
// // 	double dissOPrevious(int si);
// // 	double params(std::string name, int site);
// // private:
// // 	std::vector<double> _dissO1;
// // 	std::vector<double> _dissO2;
// // 	int _nTimeSteps;
// // 	Watershed _watershed;
// // }

// // class Site {
// // public:
// // 	std::vector<int> upstream();
// // 	double lateralWeight();
// // 	double lateralDissO();
// // 	double weight();
// // 	double area();
// // 	double dx();
// // }

// // class Watershed {
// // public:
// // 	const Site & at(int i);
// // 	int size();
// // }


// void NSModel::run() {
// 	for(int ti = 0; ti < _nTimeSteps; ++ti) {
// 		int realTime = getTime(ti);
// 		for(int si = 0; si < _watershed.size(); ++si) {
// 			double waterTemp = _watershed.waterTemperature(si, realTime);
// 			double pressure = _watershed.pressure(si, realTime);


// 			// this is done for R entry point, not for c++
// 			const Site & site = _watershed.at(si);
// 			double inputDissOMass;
// 			for(int upSiteIndex : site.upstream())
// 				inputDissOMass += _watershed.at(upSiteIndex).weight() * dissOPrevious(upSiteIndex);
// 			inputDissOMass += site.lateralWeight() * site.lateralDissO();
			
// 		}
// 	}
// }



// // stop here
// // this is enough code that it already requires testing
// // the problem of course is that classes in c++ don't place nicely with R
// // solution is to implement things in a way that R can interface, even if not as typical of C++
// // need multiple entry points to facilitate testing from R

// // at very least, need to simplify functions to take simple types when possible
// // so next step is a function that computes dDO/dt within a site; can easily test that code from R