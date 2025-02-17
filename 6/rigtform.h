#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS380_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) {
    t_ = t;
    r_ = r;
  }

  explicit RigTForm(const Cvec3& t) {
    t_ = t;
    r_ = Quat();
  }

  explicit RigTForm(const Quat& r) {
    t_ = Cvec3(0, 0, 0);
    r_ = r;
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const {
    // TODO
    Cvec4 c;
    c = r_ * a + Cvec4(t_, 0);
    return c;
  }

  RigTForm operator * (const RigTForm& a) const {
    // TODO
    RigTForm rt;
    rt.setTranslation(t_ + Cvec3(r_ * Cvec4(a.t_, 1))); // = t_ + r_ * a.t_;
    rt.setRotation(r_ * a.r_);
    return rt;
  }
};

inline RigTForm inv(const RigTForm& tform) {
  RigTForm rt;
  rt.setTranslation(Cvec3(inv(tform.getRotation()) * Cvec4(tform.getTranslation(), 0)) * -1);
  rt.setRotation(inv(tform.getRotation()));
  return rt;
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline RigTForm makeMixedFrame(const RigTForm&objectRbt, const RigTForm&eyeRbt) {
  return transFact(objectRbt) * linFact(eyeRbt);
}

inline RigTForm doMtoOwrtA(const RigTForm&M, const RigTForm&O, const RigTForm&A) {
  return A * M * inv(A) * O;
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
  Matrix4 T = Matrix4::makeTranslation(tform.getTranslation());
  Matrix4 R = quatToMatrix(tform.getRotation());

  return T * R;
}

#endif
