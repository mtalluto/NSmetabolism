#ifndef params_h
#define params_h

#include <memory>
#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>

namespace NSM {
	// handy constants
	const unsigned MIN_PER_DAY = (24 * 60);

	/**
	 * @brief      Class for holding parameters for the metabolism model 
	*/
	class Params {
	public:
		double lP1;
		double lP2;
		double er24_20;
		double k600;

		Params(const Params &p);
		Params(double lp1, double lp2, double er, double k);
		Params(const Rcpp::NumericVector &p);
	};



	/**
	 * @brief      Encapsulated data for a segment of a river
	*/
	class RData {
	public:
		double depth; // meters
		double initial_do; // g/meter^3
		double discharge; // m^3/sec
		double width; // meters
		double length; // meters
		double lateral_do; // g/meter^3
		int id;

		double cs_area() const {return depth * width;}

		RData(double z, double ido, double q, double w, double dx, double ldo, int id) :
			depth(z), initial_do(ido), discharge(q), width(w), length(dx), lateral_do(ldo),
			id(id) {}
		RData(double z, double ido) : RData(z, ido, 0, 0, 0, 0, 0) {}
	};


	/**
	  For time-series datasets
	*/

	class Timeseries {
	protected:
		const static std::set<std::string> _allowed_variables;
		std::unordered_map<std::string, std::vector<double> > _data;
		const unsigned _dt; // size of time step, in minutes

	public:
		void insert(std::string var, std::vector<double> data);
		void insert(std::string var, const Rcpp::NumericVector &data);
		double at(std::string var, size_t time) const;
		double & at(std::string var, size_t time);
		std::vector<double> at(std::string var) const;
		unsigned dt() const {return _dt;}
		size_t nt() const {return _data.cbegin()->second.size();}
		bool empty() const noexcept {return _data.empty();}

		Timeseries(unsigned dt);
		Timeseries(unsigned dt, std::string var, std::vector<double> data);
		Timeseries(unsigned dt, std::string var, const Rcpp::NumericVector &data);
	};

	class DailyTimeseries: public Timeseries {
	public:
		double & at_min(std::string var, size_t time_min);
		DailyTimeseries() : Timeseries(MIN_PER_DAY) {}
	};

	std::vector<std::shared_ptr<Params> > param_from_r(const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
		const Rcpp::NumericVector &k600);

	std::vector<std::shared_ptr<RData> > data_from_r(const Rcpp::DataFrame &df);
}

#endif
