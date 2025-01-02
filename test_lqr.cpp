/**
 * @file test_lqr.cpp
 * @brief LQR控制器测试程序
 */

#include "lqr_controller.hpp"
#include <vector>
#include <iostream>

int main() {
    // 系统参数定义
    linalg::Matrix<double> A({{1.0}});
    linalg::Matrix<double> B({{1.0}});
    linalg::Matrix<double> Q({{1.0}});
    linalg::Matrix<double> R({{1.0}});
    linalg::Matrix<double> S({{1.0}});

    // 创建LQR控制器
    control::LQRController<double> lqr(A, B, Q, R, S);

    // 初始状态
    linalg::Matrix<double> x({{1.0}});
    
    // 存储历史数据
    std::vector<double> x_history;
    std::vector<double> u_history;
    
    // 仿真步数
    const int k_steps = 20;
    
    x_history.push_back(x(0,0));
    
    // 仿真循环
    for(int k = 0; k < k_steps; ++k) {
        // 计算控制输入
        linalg::Matrix<double> u = lqr.compute_input(x);
        u_history.push_back(u(0,0));
        
        // 更新系统状态
        x = A * x + B * u;
        x_history.push_back(x(0,0));
    }

    // 输出结果
    std::cout << "状态历史：\n";
    for(size_t i = 0; i < x_history.size(); ++i) {
        std::cout << "x[" << i << "] = " << x_history[i] << "\n";
    }
    
    std::cout << "\n控制输入历史：\n";
    for(size_t i = 0; i < u_history.size(); ++i) {
        std::cout << "u[" << i << "] = " << u_history[i] << "\n";
    }

    return 0;
} 