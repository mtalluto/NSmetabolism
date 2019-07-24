functions{





/**
  * Compute transport component
  * @param inputDO dissolved oxygen concentration of input water
  * @param outputDO dissolved oxygen concentration of focal pixel
  * @param Q discharge
  * @param area cross sectional area
  * @param dx length of pixel
*/
real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx) {
	real inputMass = Q * inputDO;
	real outputMass = Q * outputDO;
	return (-1/area) * (outputMass - inputMass)/dx;
}






// }
