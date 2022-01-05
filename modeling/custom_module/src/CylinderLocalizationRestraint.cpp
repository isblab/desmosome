#include <IMP/core/XYZ.h>
#include <IMP/desmosome/CylinderLocalizationRestraint.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/algebra/Rotation3D.h>
#include <IMP/algebra/VectorD.h>
#include <math.h>

IMPDESMOSOME_BEGIN_NAMESPACE

/* This constructor initializes the class in mode 0, the general setup. Requires the coordinates of 2 points
through which the axis of the cylinder passes through, radius, kappa and the particle list */
CylinderLocalizationRestraint::CylinderLocalizationRestraint(IMP::ParticlesTemp plist, double kappa,
double axis_x1, double axis_y1, double axis_z1, double axis_x2, double axis_y2, double axis_z2, double r) :
        Restraint(plist[0]->get_model(), "CylinderLocalizationRestraint %1%"),
        plist_(plist),
        kappa_(kappa),
        axis_x1_(axis_x1),
        axis_y1_(axis_y1),
        axis_z1_(axis_z1),
        axis_x2_(axis_x2),
        axis_y2_(axis_y2),
        axis_z2_(axis_z2),
        rSquared_(r * r) {
            mode_ = 0;
            axis_x2_ -= axis_x1_;  // To find a vector along the cylinder-axis from the two given points
            axis_y2_ -= axis_y1_;
            axis_z2_ -= axis_z1_;
            // Find the magnitude of the vector along the cylinder-axis to later make it a unit vector
            double magnitude = sqrt(axis_x2_ * axis_x2_ + axis_y2_ * axis_y2_ + axis_z2_ * axis_z2_);
            // Find the angle of the cylinder-axis with the z-axis
            double theta = acos(axis_z2_ / magnitude);
            double x = sin(theta / 2);
            // Rotation is created from a Quaternion such that the rotation rotates the cylinder-axis parallel to z
            IMP::algebra::Rotation3D rot = IMP::algebra::Rotation3D(cos(theta / 2), x * axis_y2_ / magnitude / sin(theta), -x * axis_x2_ / magnitude / sin(theta), 0);
            IMP::algebra::Vector3D vec = IMP::algebra::Vector3D(0, 0, 0);
            transformation_ = IMP::algebra::Transformation3D(rot, vec);
            // Identify one point on the transformed cylinder-axis to act as the "center" -> consistent with mode 1
            vec = transformation_.get_transformed(IMP::algebra::Vector3D(axis_x1_, axis_y1_, axis_z1_));
            center1_ = vec[0];
            center2_ = vec[1];
        }

/* This constructor initializes the class in mode 1. Requires the axis (0, 1, 2 for x, y, z) to which the 
axis of the cylinder is parallel, and the coordinates along the other two axes for any point
on the actual axis (to identify the shift of the axis), radius, kappa and the particle list */
CylinderLocalizationRestraint::CylinderLocalizationRestraint(IMP::ParticlesTemp plist, double kappa,
int axis, double center1, double center2, double r) :
        Restraint(plist[0]->get_model(), "CylinderLocalizationRestraint %1%"),
        plist_(plist),
        kappa_(kappa),
        axis_(axis),
        center1_(center1),
        center2_(center2),
        rSquared_(r * r){
            mode_ = 1;  // Given the coordinate axis, identify the coordinate numbers (0, 1, 2) of the other two axes
            if (axis_ == 0){
                coord1_ = 1;
                coord2_ = 2;
            }
            else if (axis_ == 1){
                coord1_ = 0;
                coord2_ = 2;
            }
            else {
                coord1_ = 0;
                coord2_ = 1;
            }
        }

// For mode 1: Calculate the score for a single particle
double CylinderLocalizationRestraint::getRadialDistanceSpecial(IMP::Particle* p) const {
    double val1 = IMP::core::XYZ(p).get_coordinate(coord1_);  // The coordinates of the particle
    double val2 = IMP::core::XYZ(p).get_coordinate(coord2_);
    val1 -= center1_;  // Calculate the coordinate-wise distance from the center
    val2 -= center2_;
    double radial = (val1 * val1) + (val2 * val2);  // radial = (euclidean distance from center) ^ 2
    if (radial > rSquared_){  // Continue only if the distance exceeds the radius
        double deviation = rSquared_ + radial - 2 * sqrt(rSquared_ * radial);
        // deviation = (total distance from center - radius) ^ 2
        return fabs(kappa_ * deviation);
    }
    else {
        return 0;
    }
}

// For mode 2: Calculate the score for a single particle
double CylinderLocalizationRestraint::getRadialDistanceGeneral(IMP::Particle* p) const {
    double valx = IMP::core::XYZ(p).get_coordinate(0);  // The coordinates of the particle
    double valy = IMP::core::XYZ(p).get_coordinate(1);
    double valz = IMP::core::XYZ(p).get_coordinate(2);
    // Transform all the points such that the axis is now aligned to z-axis
    IMP::algebra::Vector3D pvec = transformation_.get_transformed(IMP::algebra::Vector3D(valx, valy, valz));
    double val1 = pvec[0] - center1_;  // Calculate the coordinate-wise distance from the center
    double val2 = pvec[1] - center2_;  // Only done for x and y coordinates
    double radial = (val1 * val1) + (val2 * val2);  // radial = (euclidean distance from center) ^ 2
    if (radial > rSquared_){  // Continue only if the distance exceeds the radius
        double deviation = rSquared_ + radial - 2 * sqrt(rSquared_ * radial);
        // deviation = (total distance from center - radius) ^ 2
        return fabs(kappa_ * deviation);
    }
    else {
        return 0;
    }
}

// Sums the individual particle score and returns the total score
double CylinderLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    if (mode_ == 0){
        for (unsigned int i=0; i < plist_.size(); i++){
            score += getRadialDistanceGeneral(plist_[i]);
        }
    }
    else {
        for (unsigned int i=0; i < plist_.size(); i++){
            score += getRadialDistanceSpecial(plist_[i]);
        }
    }
    if (accum){};
    return score;
}

IMP::ModelObjectsTemp CylinderLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPDESMOSOME_END_NAMESPACE