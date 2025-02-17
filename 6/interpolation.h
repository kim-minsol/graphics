#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "rigtform.h"

namespace Interpolation {
	inline Cvec3 lerp(const Cvec3& c0, const Cvec3& c1, const double& alpha) {
		return Cvec3(c0 * (1 - alpha) + c1 * alpha);	
	}

	// Spherical interpolation of two quaternions
	// TODO: compute the spherical interpolation of two Quat inputs and return the result
	inline Quat slerp(const Quat& q0, const Quat& q1, const double& alpha) {
		Quat result = Quat();
		Quat temp = q1 * inv(q0);
		if (temp[0] < 0) {
			temp *= -1;
		} 
		double p = atan2(std::sqrt(std::pow(temp[1], 2) + std::pow(temp[2], 2) + std::pow(temp[3], 2)), temp[0]);
		double ap = alpha * p;
		
		const double threshold = 0.001;
		if (p > CS380_PI - threshold) { // Large angle
			// Use lerp instead
			result = normalize(q0 * (1.0 - alpha) + q1 * alpha );
		} else if (p < threshold) { // Small angle
			// Use lerp instead
			result = normalize(q0 * (1.0 - alpha) + q1 * alpha);
		} else {
			result[0] = cos(ap);
			result[1] = temp[1] / sin(p) * sin(ap);
			result[2] = temp[2] / sin(p) * sin(ap);
			result[3] = temp[3] / sin(p) * sin(ap);
			result = result * q0;
		}
        return result;
	}

	// Linear interpolation of two RigTForms
	// TODO: compute the linear interpolation of two RigTForm inputs and return the result
	// Note: you should use the lerp and slerp functions you implemented above
	inline RigTForm Linear(const RigTForm& rbt0, const RigTForm& rbt1, const double& alpha) {
		RigTForm interpRbt = RigTForm();
		interpRbt.setTranslation(lerp(rbt0.getTranslation(), rbt1.getTranslation(), alpha));
		interpRbt.setRotation(slerp(rbt0.getRotation(), rbt1.getRotation(), alpha));
		return interpRbt;	
	}

	// Catmull-Rom interpolation of two RigTForms
	// TODO: compute the Catmull-Rom interpolation of four RigTForm inputs and return the result
	// Note: To Catmull-Rom interpolate two RigTFrom rbt0 and rbt1, we need 4 Keyframes keyframe rbt_1, rbt0, rbt1, rbt2.
	// 		 You can use the lerp and slerp functions you implemented above.
	inline RigTForm CatmullRom(const RigTForm& rbt_1, const RigTForm& rbt0, const RigTForm& rbt1, const RigTForm& rbt2, const double& alpha) {
		// Translation interp
		Cvec3 ci_1 = rbt_1.getTranslation();
		Cvec3 ci = rbt0.getTranslation();
		Cvec3 ci1 = rbt1.getTranslation();
		Cvec3 ci2 = rbt2.getTranslation();
		
		Cvec3 di = (ci1 - ci_1) / 6 + ci;
		Cvec3 ei = (ci2 - ci) / 6 * -1. + ci1;

		Cvec3 p01 = lerp(ci, di, alpha);
		Cvec3 p12 = lerp(di, ei, alpha);
		Cvec3 p23 = lerp(ei, ci1, alpha);
		Cvec3 p012 = lerp(p01, p12, alpha);
		Cvec3 p123 = lerp(p12, p23, alpha);
		Cvec3 p = lerp(p012, p123, alpha);

		// Rotation interp
		Quat qi_1 = rbt_1.getRotation();
		Quat qi = rbt0.getRotation();
		Quat qi1 = rbt1.getRotation();
		Quat qi2 = rbt2.getRotation();

		Quat qdi = pow(qi1*inv(qi_1), (1/6))*qi;
		Quat qei = pow(qi2*inv(qi), (-1/6))*qi1;
		
		Quat q01 = slerp(qi, qdi, alpha);
		Quat q12 = slerp(qdi, qei, alpha);
		Quat q23 = slerp(qei, qi1, alpha);
		Quat q012 = slerp(q01, q12, alpha);
		Quat q123 = slerp(q12, q23, alpha);
		Quat q = slerp(q012, q123, alpha); 

		
		return RigTForm(p, q);	// Replace this value with your own code
	}
}
#endif