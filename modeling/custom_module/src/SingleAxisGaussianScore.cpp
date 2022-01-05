/* Restraint to calculate a Gaussian-score based restraint using a single axis coordinate of the particles */


#include <IMP/desmosome/SingleAxisGaussianScore.h>  //To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/core/XYZ.h>

IMPDESMOSOME_BEGIN_NAMESPACE

SingleAxisGaussianScore::SingleAxisGaussianScore(int axis, double axisMean, double sigma) :
		sigma_(sigma), axis_(axis), axisMean_(axisMean) {}
		

double SingleAxisGaussianScore::evaluate_index(IMP::Model* m, IMP::ParticleIndex pi, IMP::DerivativeAccumulator* accum) const {
    double coord = core::XYZ(m, pi).get_coordinate(axis_);
    double score = 1. / (sqrt(2. * IMP::PI) * sigma_) * exp(-pow((coord - axisMean_), 2.) / (2. * sigma_ * sigma_));
    if (accum){};
    return -log(score);
}
//IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
IMP::ModelObjectsTemp SingleAxisGaussianScore::do_get_inputs(IMP::Model* m, const IMP::ParticleIndexes &pis) const {
    return IMP::get_particles(m, pis);
}


IMPDESMOSOME_END_NAMESPACE

