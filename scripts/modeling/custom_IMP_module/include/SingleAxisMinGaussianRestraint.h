/* Restraint to calculate a Gaussian-score based restraint using a single axis coordinate of the particles */

#ifndef IMPDESMOSOME_SINGLE_AXIS_MIN_GAUSSIAN_RESTRAINT_H  // Header guard
#define IMPDESMOSOME_SINGLE_AXIS_MIN_GAUSSIAN_RESTRAINT_H

#include <IMP/desmosome/desmosome_config.h>  //To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>

IMPDESMOSOME_BEGIN_NAMESPACE

class IMPDESMOSOMEEXPORT SingleAxisMinGaussianRestraint : public IMP::Restraint {
	double sigma_;  //the sigma for the gaussian
	int axis_;  //axis identifier (0, 1, 2)
	double axisMean_;  //the mean of the gaussian along the axis
	IMP::ParticlesTemp plist_;  //all the particles for which the minimum score is calculated
	
	public:
		SingleAxisMinGaussianRestraint(IMP::ParticlesTemp plist, int axis, double axisMean, double sigma);
		
		// unprotected_evaluate calculates the score
		// do_get_inputs returns the particles for which the score was calculated
		virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
		//IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
		virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
		IMP_OBJECT_METHODS(SingleAxisMinGaussianRestraint);  //add the usual IMP object methods
	
	private:
		double getMinCoord() const;
};

IMPDESMOSOME_END_NAMESPACE

#endif