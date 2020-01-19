#ifndef params_h
#define params_h

#include <memory>

/**
 * @brief      Class for holding parameters for the metabolism model 
*/
namespace NSM {
	class Params {
	public:
		double lP1;
		double lP2;
		double er24_20;
		double k600;
	};

	typedef std::shared_ptr<Params> param_ptr;

}

#endif
