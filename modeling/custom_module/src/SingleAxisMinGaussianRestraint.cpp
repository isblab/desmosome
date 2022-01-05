#include <IMP/core/XYZ.h>
#include <IMP/desmosome/SingleAxisMinGaussianRestraint.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <limits>

IMPDESMOSOME_BEGIN_NAMESPACE

SingleAxisMinGaussianRestraint::SingleAxisMinGaussianRestraint(IMP::ParticlesTemp plist, int axis, double axisMean, double sigma) :
		Restraint(plist[0]->get_model(), "SingleAxisMinGaussianRestraint %1%"),
		sigma_(sigma),
		axis_(axis),
		axisMean_(axisMean),
		plist_(plist) {}

double SingleAxisMinGaussianRestraint::getMinCoord() const {
	double currentMin =  std::numeric_limits<double>::infinity();
	double currVal = 0;
	for (long unsigned int i = 0; i < plist_.size(); i++){  //find the minimum distance of the particle from the mean
		//take the coordinate using the decorator XYZ
		currVal = core::XYZ(plist_[i]).get_coordinate(axis_);
		if (fabs(currVal - axisMean_) <= fabs(currentMin - axisMean_)){
			currentMin = currVal;
		}
	}
	return currentMin;
}

double SingleAxisMinGaussianRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
	double minCoord = getMinCoord();
	double score = 1. / (sqrt(2. * IMP::PI) * sigma_) * exp(-pow((minCoord - axisMean_), 2.) / (2. * sigma_ * sigma_));
	if (accum){};
	return -log(score);
}

IMP::ModelObjectsTemp SingleAxisMinGaussianRestraint::do_get_inputs() const {
	return plist_;
}


IMPDESMOSOME_END_NAMESPACE