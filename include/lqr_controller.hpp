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
    const linalg::Matrix<T>& get_gain() const { 
        return F; 
    }
};

template<typename T>
linalg::Matrix<T> LQRController<T>::calculate_gain() {
    // 获取系统维度
    size_t n = A.get_rows();    // 状态维度
    size_t p = B.get_cols();    // 输入维度
    
    // 初始化P0为终值权重矩阵S
    linalg::Matrix<T> P_k_min_1 = S;
    
    // 定义迭代参数
    const size_t max_iter = 200;
    const T tol = T(1e-3);      // 收敛阈值
    T diff = std::numeric_limits<T>::infinity();
    
    // 初始化增益矩阵
    linalg::Matrix<T> F_N_min_k(p, n);     // 当前步增益
    linalg::Matrix<T> F_N_min_k_pre(p, n); // 上一步增益
    
    // 迭代计算
    size_t k = 1;
    while(diff > tol) {
        // 保存上一步的增益矩阵
        F_N_min_k_pre = F_N_min_k;
        
        // 计算F[N-k] = (R + B'P[k-1]B)^(-1) * B'P[k-1]A
        linalg::Matrix<T> Bt = B.transpose();
        linalg::Matrix<T> temp1 = Bt * P_k_min_1;
        linalg::Matrix<T> temp2 = temp1 * B;
        F_N_min_k = (temp2 + R).inverse() * temp1 * A;
        
        // 计算P[k] = (A-BF)'P[k-1](A-BF) + F'RF + Q
        linalg::Matrix<T> BF = B * F_N_min_k;
        linalg::Matrix<T> A_BF = A - BF;
        linalg::Matrix<T> A_BF_t = A_BF.transpose();
        
        linalg::Matrix<T> term1 = A_BF_t * P_k_min_1 * A_BF;
        linalg::Matrix<T> F_t = F_N_min_k.transpose();
        linalg::Matrix<T> term2 = F_t * R * F_N_min_k;
        linalg::Matrix<T> P_k = term1 + term2 + Q;
        
        // 更新P[k-1]
        P_k_min_1 = P_k;
        
        // 计算收敛误差
        diff = (F_N_min_k - F_N_min_k_pre).max().abs().scalar();
        // 迭代计数
        ++k;
        if(k > max_iter) {
            throw std::runtime_error("Maximum Number of Iterations Exceeded");
        }
    }
    
    // 返回最终的增益矩阵
    return F_N_min_k;
}

} // namespace control