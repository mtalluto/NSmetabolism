#include "../inst/include/reach.h"
#include "../inst/include/params.h"
#include "../inst/include/pixel.h"
#include <Rcpp.h>

using std::vector;
using std::shared_ptr;
using std::make_shared;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;


std::shared_ptr<NSM::Pixel> NSM::Reach::bottom() {
	if(!_topology_built)
		std::runtime_error("Tried to access a non-built _topology");
	return _bottom;
}

// constructor for the case where light, pressure, and temperature are constant across the reach
// NSM::Reach::Reach(const NumericVector &light, const NumericVector &temperature, 
// 		const NumericVector &pressure, const NSM::Params & pars) : Reach(pars)
// {
// 	int nt = light.size();
// 	if(temperature.size() != nt || pressure.size() != nt)
// 		throw std::length_error("Light, temperature, and pressure must have the same length");

// 	for(int i = 0; i < nt; ++i) {
// 		_light->push_back(light.at(i));
// 		_temperature->push_back(temperature.at(i));
// 		_pressure->push_back(pressure.at(i));
// 	}
// }


// constructor for the case where light, pressure, and temperature vary by pixel
NSM::Reach::Reach(shared_ptr<NSM::Params> pars, const Rcpp::DataFrame &pixData, 
			vector<shared_ptr<NSM::Timeseries> > ts_data, const Rcpp::IntegerMatrix &topology) : 
		Reach(pars)
{
	if(ts_data.size() != pixData.nrow())
		throw std::length_error("Timeseries must have one entry per pixel");

	// First pass, set up pixels
	auto data = data_from_r(pixData);
	for(int p = 0; p < pixData.nrow(); ++p)
		_pixels[data.at(p)->id] = make_shared<NSM::Pixel> (_pars, data.at(p), ts_data.at(p));

	// second pass, set topology
	set_topology(topology);
}


NSM::Reach::Reach(shared_ptr<Params> p) : _pars(p), _topology_built(false)
{ }


/*
	Topology matrix T is always two columns and indicates the coordinates of the adjacency matrix
	M that are non zero; so a row of T == c(i,j indicates that M[i,j] is nonzero
	This means that pixel i receives water from pixel j
	This also means that:
		-confluences will have multiple (generally 2) entries in the first column of T
		-headwaters have no entries in column 1
		-no pixel should appear more than once in column 2
		-the outlet will have no entry in column 2
	Additionally, because reaches can be connected to other reaches, we have to check for
	pixels that have an entry in the topology matrix but no entry in the pixel structure; in 
	this case we search in the neighboring reaches
*/
void NSM::Reach::set_topology(const Rcpp::IntegerMatrix &T) {
	for(size_t i = 0; i < T.nrow(); ++i) {
		std::unordered_map<int, shared_ptr<Pixel>>::iterator downstream = _pixels.find(T(i, 0));
		std::unordered_map<int, shared_ptr<Pixel>>::iterator upstream = _pixels.find(T(i, 1));

		// 4 possibilities, neither found, us found but ds not, ds found but us not, both found
		if(downstream == _pixels.end()) {
			if(upstream == _pixels.end()) {
				// not found, it's irrelevant, only here to be explicit
				continue;
			} else {
				// downstream in another reach, but upstream is present; this must be the bottom
				// of the reach
				if(_bottom != nullptr && _bottom != upstream->second)
					throw std::runtime_error("Invalid topology");
				_bottom = upstream->second;
			}
		} else {
			// downstream is in this reach
			// if upstream not found, we must be at the top, need to find the neighbor in
			// the next reach up
			if(upstream == _pixels.end()) {
				if(_top != nullptr)
					throw std::runtime_error("Invalid topology");
				_top = downstream->second;
				for(auto it : _upstream) {
					if(it->bottom()->id() == T(i, 1))
						downstream->second->add_upstream(it->bottom());
				}
			} else {
				// both pixels in this reach
				downstream->second->add_upstream(upstream->second);
			}
		}
	}

	// note that if the top of the reach is a headwater, it won't have been detected so
	// needs to be found manually; same for the outlet; outlet will have its upstream set,
	// but it won't have been id'ed as the bottom yet
	// outlet is the pixel that is not upstream of anyone
	// headwaters are not downstream of anyone
	auto p = _pixels.begin();
	vector<int> uspix = Rcpp::as<vector<int> >(Rcpp::IntegerVector(T.row(1)));
	while((_bottom == nullptr || _top == nullptr) && p != _pixels.end()) {
		if(p->second->upstream().empty() && _top == nullptr)
			_top = p->second;
		vector<int>::iterator it = std::find(uspix.begin(), uspix.end(), p->first);
		if(it == uspix.end())
			_bottom = p->second;
		++p;
	}
	_topology_built = true;
}
