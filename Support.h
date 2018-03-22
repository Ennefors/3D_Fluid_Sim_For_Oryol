//
// Created by bill on 19/01/18.
//

#ifndef ORYOL_SUPPORT_H
#define ORYOL_SUPPORT_H
#pragma once

#include <glm/vec3.hpp>

#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <Eigen/Sparse>




namespace support
{


    inline double rk4(double(*f)(double, double), double dx, double x, double y)
    {
        double  k1 = dx * f(x, y),
                k2 = dx * f(x + dx / 2, y + k1 / 2),
                k3 = dx * f(x + dx / 2, y + k2 / 2),
                k4 = dx * f(x + dx, y + k3);
        return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    //Calculates the gradient of the given vector
    inline glm::vec3 calc_gradient(float* g(float, float , float), glm::vec3 v)
    {
        return glm::vec3(g(v.x, v.y, v.z) - g(v.x-1, v.y, v.z),
                         g(v.x, v.y, v.z) - g(v.x, v.y-1, v.z),
                         g(v.x, v.y, v.z) - g(v.x, v.y, v.z-1));
    }

    inline float calc_laplace(float g(float x, float y, float z, int ind), glm::vec3 v, int ind)
    {

        float vx = g(v.x+1, v.y, v.z, ind) + g(v.x-1, v.y, v.z,ind);
        float vy = g(v.x, v.y+1, v.z, ind) + g(v.x, v.y-1, v.z, ind);
        float vz = g(v.x, v.y, v.z+1, ind) + g(v.x, v.y, v.z-1, ind);
        float v0 = g(v.x, v.y, v.z, ind)*6;

        return ((vx) + (vy) + (vz) - v0);
    }

    inline float calc_CpartDeriv(float* g(float, float, float), glm::vec3 v)
    {
        return g(v.x, v.y+1, v.z) - g(v.x, v.y-1, v.z);
    }

    inline float calc_FpartDeriv(float* g(float, float, float), glm::vec3 v)
    {
        return g(v.x, v.y+1, v.z) - g(v.x, v.y, v.z);
    }

    inline float calc_BpartDeriv(float* g(float, float, float), glm::vec3 v)
    {
        return g(v.x, v.y, v.z) - g(v.x, v.y-1, v.z);
    }

    inline float calc_divVector(glm::vec3 v)
    {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;

        float vxf = v.x+1;
        float vyf = v.y+1;
        float vzf = v.z+1;

        return ((vxf - vx) +
                (vyf - vy) +
                (vzf - vz));
    }


    template<class T, class Compare>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi, Compare comp )
    {
        return (comp(v, lo) ? lo : comp(hi, v) ? hi : v);
    }

    template<class T>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi )
    {
        return clamp( v, lo, hi, std::less<float>() );
    }


}
#endif //ORYOL_SUPPORT_H
