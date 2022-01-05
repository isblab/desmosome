/* Restraint to calculate a Harmonic restraint for particles if they exit a cylindrical boundary
about an arbitrary axis defined by passing through two points x1, y1, z1 and x2, y2, z2. If the axis
aligns with one of the co-ordinate axes, can alternatively initialize the class to save some computation */

#ifndef IMPDESMOSOME_CYLINDER_LOCALIZATION_RESTRAINT_H  // Header guard
#define IMPDESMOSOME_CYLINDER_LOCALIZATION_RESTRAINT_H

#include <IMP/desmosome/desmosome_config.h>  // To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>
#include <IMP/algebra/Transformation3D.h>

IMPDESMOSOME_BEGIN_NAMESPACE

class IMPDESMOSOMEEXPORT CylinderLocalizationRestraint : public IMP::Restraint {
    IMP::ParticlesTemp plist_; // All the particles to which to apply this to (individual scores are summed)
    double kappa_;  // The kappa for the harmonic
    double axis_x1_;  // From the equation of a line t(x2-x1, y2-y1, z2-z1)
    double axis_y1_;  // The axis passes through the two points x1, y1, z1 and x2, y2, z2
    double axis_z1_;
    double axis_x2_;
    double axis_y2_;
    double axis_z2_;
    int axis_;  // Can alternatively just specify the axis if the cylinder is parallel to one of the coordinate axes
    double center1_;  // If axis_ is given, give the coordinates of one point on the axis (otherwise internal only)
    double center2_;  // The order of the coordinates is x -> y -> z with the coordinate along axis_ skipped
    double rSquared_;  // Stores the square of the radius of the cylinder
    int mode_;  // 0 for the general initialization, 1 for the initialization with axis-parallel to x, y or z
    int coord1_;  // (internal only) To store the axis indices (0, 1, 2) for mode 1
    int coord2_;
    IMP::algebra::Transformation3D transformation_;  
    // (internal only) For mode 0, a transformation that rotates the axis of the cylinder to be aligned to z-axis

    public:
        CylinderLocalizationRestraint(IMP::ParticlesTemp plist, double kappa, double axis_x1,
        double axis_y1, double axis_z1, double axis_x2, double axis_y2, double axis_z2, double r);
        CylinderLocalizationRestraint(IMP::ParticlesTemp plist, double kappa, int axis, double center1,
        double center2, double r);

        // unprotected_evaluate calculates the score
        // do_get_inputs returns the particles for which the score was calculated
        virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
        //IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
        virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
        IMP_OBJECT_METHODS(CylinderLocalizationRestraint);  //add the usual IMP object methods

    private:
        double getRadialDistanceGeneral(IMP::Particle* p) const;
        double getRadialDistanceSpecial(IMP::Particle* p) const;
};

IMPDESMOSOME_END_NAMESPACE

#endif