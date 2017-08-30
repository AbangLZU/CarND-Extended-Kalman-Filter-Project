#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
    TODO:
      * predict the state
    */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
      * update the state by using Kalman Filter equations
    */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;

    KF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
      * update the state by using Extended Kalman Filter equations
    */

    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];

    float rho = sqrt(px * px + py * py);
    float phi = atan2(py, px);
    float rho_dot;
    if(fabs(rho) < 0.00001){
        rho_dot = 0.0;
    } else{
        rho_dot =(px * vx + py * vy)/rho;
    }

    VectorXd hx = VectorXd(3);
    hx << rho, phi, rho_dot;
    VectorXd y = z - hx;
    KF(y);
}

void KalmanFilter::KF(Eigen::VectorXd &y) {
    // normalize the y to -pi~pi
    double const PI=3.14159265;
    while(y[1] > PI){
        y[1] = y[1] - 2 * PI;
    }
    while(y[1] < -PI){
        y[1] = y[1] + 2 * PI;
    }

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
