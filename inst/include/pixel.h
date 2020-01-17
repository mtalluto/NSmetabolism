#ifndef pixel_h
#define pixel_h

#include <memory>
#include <vector>
#include <Rcpp.h>
#include "../inst/include/params.h"

namespace NSM {
	// std::vector<Pixel> dfToPixel(const Rcpp::DataFrame &pixDf);

	// class Pixel {
	// 	static int nt;  // number of time steps, shared by all Pixels
	// 	static int dt;	// size of time steps, in seconds

	// 	std::shared_ptr<Params> pars;
	// 	int id;
	// 	double discharge;
	// 	// double velocity;
	// 	// double depth;
	// 	// double width;
	// 	// double length;

	// 	// eventually provide a friend class LightArray that does the heavy lifting;
	// 	// accessed with light(id, time); thus can abstract all the spatial crap
	// 	double light [nt]; 

	// 	// eventually
	// 	// std::shared_ptr<Pixel> downstream;
	// 	// std::vector<std::shared_ptr<Pixel> > upstream;
	// public:
	// 	double gpp(int t);
	// };
}

#endif
