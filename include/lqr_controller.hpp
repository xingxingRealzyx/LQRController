/**
 * @file lqr_controller.hpp
 * @brief LQR控制器实现
 * @author xingxing
 * @date 2025
 */

#pragma once
#include "linalg.hpp"
#include <cmath>

namespace control {

/**
 * @brief LQR控制器类
 * @tparam T 数值类型
 */
template<typename T>
class LQRController {
private:
    linalg::Matrix<T> A;  ///< 系统矩阵
    linalg::Matrix<T> B;  ///< 输入矩阵
    linalg::Matrix<T> Q;  ///< 状态权重矩阵
    linalg::Matrix<T> R;  ///< 输入权重矩阵
    linalg::Matrix<T> S;  ///< 终值权重矩阵
    linalg::Matrix<T> F;  ///< 反馈增益矩阵

    /**
     * @brief 计算LQR反馈增益矩阵
     * @return 计算得到的反馈增益矩阵
     */
    linalg::Matrix<T> calculate_gain();

public:
    /**
     * @brief 构造函数
     * @param A_in 系统矩阵
     * @param B_in 输入矩阵
     * @param Q_in 状态权重矩阵
     * @param R_in 输入权重矩阵
     * @param S_in 终值权重矩阵
     */
    LQRController(const linalg::Matrix<T>& A_in, 
                  const linalg::Matrix<T>& B_in,
                  const linalg::Matrix<T>& Q_in,
                  const linalg::Matrix<T>& R_in,
                  const linalg::Matrix<T>& S_in)
        : A(A_in), B(B_in), Q(Q_in), R(R_in), S(S_in), 
          F(linalg::Matrix<T>(B_in.get_cols(), A_in.get_cols())) {
        F = calculate_gain();
    }

    /**
     * @brief 计算控制输入
     * @param x 当前状态
     * @return 计算得到的控制输入
     */
    linalg::Matrix<T> compute_input(const linalg::Matrix<T>& x) {
        linalg::Matrix<T> result = F * x;
        for(size_t i = 0; i < result.get_rows(); ++i) {
            for(size_t j = 0; j < result.get_cols(); ++j) {
                result(i,j) *= T(-1);
            }
        }
        return result;
    }

    /**
     * @brief 获取反馈增益矩阵
     * @return 当前的反馈增益矩阵
     */
    const linalg::Matrix<T>& get_gain() const { return F; }
};

template<typename T>
linalg::Matrix<T> LQRController<T>::calculate_gain() {
    // 对于一维系统，直接计算稳态解
    T s = (Q(0,0) + std::sqrt(Q(0,0)*Q(0,0) + 4*Q(0,0)*R(0,0))) / 2;
    
    // 计算反馈增益
    linalg::Matrix<T> Bt = B.transpose();
    linalg::Matrix<T> sB = B;
    // 将B矩阵乘以标量s
    for(size_t i = 0; i < sB.get_rows(); ++i) {
        for(size_t j = 0; j < sB.get_cols(); ++j) {
            sB(i,j) *= s;
        }
    }
    // 计算(B^T * s * B + R)^(-1)
    linalg::Matrix<T> temp = (Bt * sB + R).inverse();
    
    // 计算最终的增益矩阵：K = (B^T * s * B + R)^(-1) * B^T * s * A
    linalg::Matrix<T> sA = A;
    // 将A矩阵乘以标量s
    for(size_t i = 0; i < sA.get_rows(); ++i) {
        for(size_t j = 0; j < sA.get_cols(); ++j) {
            sA(i,j) *= s;
        }
    }
    return temp * Bt * sA;
}

} // namespace control