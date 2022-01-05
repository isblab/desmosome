/* Restraint to calculate a Gaussian-score based restraint using a single axis coordinate of the particles */

#ifndef IMPDESMOSOME_SINGLE_AXIS_GAUSSIAN_SCORE_H  // Header guard
#define IMPDESMOSOME_SINGLE_AXIS_GAUSSIAN_SCORE_H

#include <IMP/desmosome/desmosome_config.h>  //To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>


IMPDESMOSOME_BEGIN_NAMESPACE

class IMPDESMOSOMEEXPORT SingleAxisGaussianScore : public IMP::SingletonScore {
	double sigma_;  //the sigma for the gaussian
	int axis_;  //axis identifier (0, 1, 2)
	double axisMean_;  //the mean of the gaussian along the axis
	
	public:
		SingleAxisGaussianScore(int axis, double axisMean, double sigma);
		
		// evaluate_index calculates the score
		// do_get_inputs returns the particles for which the score was calculated
		virtual double evaluate_index(IMP::Model* m, IMP::ParticleIndex p, IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
		//IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
		virtual IMP::ModelObjectsTemp do_get_inputs(IMP::Model* m, const IMP::ParticleIndexes &pis) const IMP_OVERRIDE;
		IMP_SINGLETON_SCORE_METHODS(SingleAxisGaussianScore);  // add the Singleton_score methods
		IMP_OBJECT_METHODS(SingleAxisGaussianScore);  //add the usual IMP object methods
	
	private:
};

IMPDESMOSOME_END_NAMESPACE

#endif