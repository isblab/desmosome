#include <IMP/core/XYZ.h>
#include <IMP/desmosome/AxialLocalizationRestraint.h>
#include <math.h>

IMPDESMOSOME_BEGIN_NAMESPACE

AxialLocalizationRestraint::AxialLocalizationRestraint(IMP::ParticlesTemp plist, int axis, double axisUpper,
double axisLower, double kappa) :
        Restraint(plist[0]->get_model(), "AxialLocalizationRestraint %1%"),
        plist_(plist),
        axis_(axis),
        axisUpper_(axisUpper),
        axisLower_(axisLower),
        kappa_(kappa){}

// Calculate a single particle score
double AxialLocalizationRestraint::getSingleScore(IMP::Particle* p) const {
    double val = core::XYZ(p).get_coordinate(axis_);  // Get the particle coordinate along the axis
    double score = 0;
    if (val > axisUpper_){  // If it exceeds the upper limit
        double x = (axisUpper_ - val) * (axisUpper_ - val);
        score += x * kappa_;
    }
    else if (val < axisLower_){  // If it exceeds the lower limit
        double x = (val - axisLower_) * (val - axisLower_);
        score += x * kappa_;
    }
    return score;
}

double AxialLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double finalScore = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        finalScore += getSingleScore(plist_[i]);
    }
    if (accum){};
    return finalScore;
}

IMP::ModelObjectsTemp AxialLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPDESMOSOME_END_NAMESPACE