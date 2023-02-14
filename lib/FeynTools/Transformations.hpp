#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include "Vectors.hpp"

using namespace Vectors;

FourTensor LorentzBoostX(double khi);
FourTensor LorentzBoostY(double khi);
FourTensor LorentzBoostZ(double khi);

FourTensor RotationX(double theta);
FourTensor RotationY(double theta);
FourTensor RotationZ(double theta);

/**
 * @brief Create a FourTensor that transforms into the rest frame of a gizen FourVector
 * 
 * @param p a FourVector, we are transfroming to its rest frame 
 * @return a FourTensor representing a Lorentz transformation to the rest frame of p.
 */
FourTensor TransformationTo(FourVector p);
FourTensor TransformationFrom(FourVector p);

#endif // TRANSFORMATIONS_HPP
