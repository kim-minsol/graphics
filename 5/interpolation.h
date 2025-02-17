#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "rigtform.h"

namespace Interpolation {

	// Linear interpolation of two coordinate vectors
	// TODO: compute the linear interpolation of two Cvec3 inputs and return the result
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
		// Quat q = q1 * inv(q0);
		// if (q[0]< 0) {
		// 	q = q * -1;
		// }
		// Quat qq1 = Quat();
		// double dotProduct = dot(q0, q1);
		// if (dotProduct < 0.0) {
		// 	qq1[0] = -q1[0];
		// 	qq1[1] = -q1[1];
		// 	qq1[2] = -q1[2];
		// 	qq1[3] = -q1[3];
		// 	dotProduct = -dotProduct;
		// }

		// //float p = atan2(sqrt(pow(q[1], 2) + pow(q[2], 2) + pow(q[3], 2)), q[0]);

		// double theta_0 = atan2(sqrt(pow(q[1], 2) + pow(q[2], 2) + pow(q[3], 2)), q[0]); // 초기 각도
		// double theta = theta_0 * alpha;             // 보간된 각도
		// double sin_theta = sin(theta);
		// double sin_theta_0 = sin(theta_0);

		// double s1 = cos(theta) - dotProduct * sin_theta / sin_theta_0;
		// double s2 = sin_theta / sin_theta_0;

		// return Quat(
		// 	s1 * q0[0] + s2 * qq1[0],
		// 	s1 * q0[1] + s2 * qq1[1],
		// 	s1 * q0[2] + s2 * qq1[2],
		// 	s1 * q0[3] + s2 * qq1[3]
		// );
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

}
#endif