//
//  Vectors with Minkowski space
//

#ifndef RSYMSQCD_VEC4D_HPP
#define RSYMSQCD_VEC4D_HPP

template<typename T>
class Vec4D {
public:
    // @todo:   understand why the list initilizer of vector of Vec4D objects needs an empty constructor
    Vec4D() {};
    explicit Vec4D(T t, T x, T y, T z) : t_(t), x_(x), y_(y), z_(z) {};
    T t_, x_, y_, z_;
};

template <typename T>
inline Vec4D<T> operator+(const Vec4D<T>& lhs, const Vec4D<T>& rhs) {
    return Vec4D<T> { lhs.t_ + rhs.t_, lhs.x_ + rhs.x_, lhs.y_ + rhs.y_, lhs.z_ + rhs.z_ };
}

template <typename T>
inline Vec4D<T> operator-(const Vec4D<T>& lhs, const Vec4D<T>& rhs) {
    return Vec4D<T> { lhs.t_ - rhs.t_, lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_ };
}

template <typename T>
inline T operator*(const Vec4D<T>& lhs, const Vec4D<T>& rhs) {
    return lhs.t_*rhs.t_ - lhs.x_*rhs.x_ - lhs.y_*rhs.y_ - lhs.z_*rhs.z_;
}

template <typename T>
inline Vec4D<T> operator*(const Vec4D<T>& lhs, T rhs) {
    return Vec4D<T> { lhs.t_*rhs, lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs };
}

template <typename T>
inline Vec4D<T> operator*(T lhs, const Vec4D<T>& rhs) {
    return Vec4D<T> { lhs*rhs.t_, lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_ };
}
#endif //RSYMSQCD_VEC4D_HPP
