/* Restraint to calculate a Harmonic restraint using a single axis coordinate of the particles.
 Penalizes the particle when it goes out of given upper and lower bounds */

#ifndef IMPDESMOSOME_SINGLE_AXIAL_LOCALIZATION_RESTRAINT_H  // Header guard
#define IMPDESMOSOME_SINGLE_AXIAL_LOCALIZATION_RESTRAINT_H

#include <IMP/desmosome/desmosome_config.h>  //To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>

IMPDESMOSOME_BEGIN_NAMESPACE

class IMPDESMOSOMEEXPORT AxialLocalizationRestraint : public IMP::Restraint {
    IMP::ParticlesTemp plist_;  // All the particles for which the minimum score is calculated
    int axis_;  // Axis identifier (0, 1, 2 for x, y, z)
    double axisUpper_;  // The upper limit along the axis
    double axisLower_;  // The lower limit
    double kappa_;  // Kappa for the harmonic


    public:
        AxialLocalizationRestraint(IMP::ParticlesTemp plist, int axis, double axisUpper, double axisLower, double kappa);
        
        // unprotected_evaluate calculates the score
        // do_get_inputs returns the particles for which the score was calculated
        virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
        //IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
        virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
        IMP_OBJECT_METHODS(AxialLocalizationRestraint);  //add the usual IMP object methods

    private:
        double getSingleScore(IMP::Particle* p) const;
};

IMPDESMOSOME_END_NAMESPACE

#endif