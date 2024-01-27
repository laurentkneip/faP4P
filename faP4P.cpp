#include "faP4P.hpp"
#include "threeConics.hpp"

Eigen::Matrix3d
polyjam::faP4P::composeB( std::vector< Eigen::Vector3d> & x ) {

  Eigen::Matrix3d temp1, temp2, temp3;
  temp1.col(0) = x[3]; temp1.col(1) = x[1]; temp1.col(2) = x[2];
  temp2.col(0) = x[0]; temp2.col(1) = x[3]; temp2.col(2) = x[2];
  temp3.col(0) = x[0]; temp3.col(1) = x[1]; temp3.col(2) = x[3];

  Eigen::Matrix3d B;
  B.col(0) = temp1.determinant() * x[0];
  B.col(1) = temp2.determinant() * x[1];
  B.col(2) = temp3.determinant() * x[2];

  return B;
}

Eigen::Matrix3d
polyjam::faP4P::composeG( std::vector< Eigen::Vector3d> & X ) {

  Eigen::Matrix3d G;

  G.col(0) = X[0]-X[3];
  G.col(1) = X[1]-X[3];
  G.col(2) = X[2]-X[3];

  return G;

}

Eigen::Matrix3d
polyjam::faP4P::getA( Eigen::Matrix3d & G ) {

  Eigen::Matrix3d C;
  for( int i = 0; i < 3; i++ ) {
    for( int j = 0; j < 3; j++ ) {
      Eigen::Matrix3d minor = G;
      minor.row(i).setZero();   // Set row i to zeros
      minor.col(j).setZero();   // Set column j to zeros
      minor(i, j) = 1.0;
      C(i,j) = minor.determinant();
    }
  }

  Eigen::Matrix3d Gadj = C.transpose();

  Eigen::Matrix3d A = Gadj * Gadj.transpose();
  return A;
}

void
polyjam::faP4P::createCoeffs( Eigen::Matrix3d & A, Eigen::Matrix3d & B,
                              Eigen::Matrix3d & A1, Eigen::Vector3d & b1, double & c1,
                              Eigen::Matrix3d & A2, Eigen::Vector3d & b2, double & c2,
                              Eigen::Matrix3d & A3, Eigen::Vector3d & b3, double & c3 ) {

  A1(0,0) = B(1,0)*B(0,0)*A(0,0);
  A1(0,1) = 0.5 * (B(1,1)*B(0,0)*A(0,1)+B(1,0)*B(0,1)*A(1,0));
  A1(1,1) = B(1,1)*B(0,1)*A(1,1);
  A1(0,2) = 0.5 * (B(1,2)*B(0,0)*A(0,2)+B(1,0)*B(0,2)*A(2,0));
  A1(1,2) = 0.5 * (B(1,2)*B(0,1)*A(1,2)+B(1,1)*B(0,2)*A(2,1));
  A1(2,2) = B(1,2)*B(0,2)*A(2,2);
  b1[0] = B(1,2)*B(0,0)*A(0,1)+B(1,2)*B(0,0)*A(0,0)+B(1,1)*B(0,0)*A(0,2)+B(1,1)*B(0,0)*A(0,0)+B(1,0)*B(0,2)*A(1,0)+B(1,0)*B(0,2)*A(0,0)+B(1,0)*B(0,1)*A(2,0)+B(1,0)*B(0,1)*A(0,0)+B(1,0)*B(0,0)*A(2,0)+B(1,0)*B(0,0)*A(1,0)+B(1,0)*B(0,0)*A(0,2)+B(1,0)*B(0,0)*A(0,1);
  b1[1] = B(1,2)*B(0,1)*A(1,1)+B(1,2)*B(0,1)*A(1,0)+B(1,1)*B(0,2)*A(1,1)+B(1,1)*B(0,2)*A(0,1)+B(1,1)*B(0,1)*A(2,1)+B(1,1)*B(0,1)*A(1,2)+B(1,1)*B(0,1)*A(1,0)+B(1,1)*B(0,1)*A(0,1)+B(1,1)*B(0,0)*A(2,1)+B(1,1)*B(0,0)*A(1,1)+B(1,0)*B(0,1)*A(1,2)+B(1,0)*B(0,1)*A(1,1);
  b1[2] = B(1,2)*B(0,2)*A(2,1)+B(1,2)*B(0,2)*A(2,0)+B(1,2)*B(0,2)*A(1,2)+B(1,2)*B(0,2)*A(0,2)+B(1,2)*B(0,1)*A(2,2)+B(1,2)*B(0,1)*A(0,2)+B(1,2)*B(0,0)*A(2,2)+B(1,2)*B(0,0)*A(1,2)+B(1,1)*B(0,2)*A(2,2)+B(1,1)*B(0,2)*A(2,0)+B(1,0)*B(0,2)*A(2,2)+B(1,0)*B(0,2)*A(2,1);
  c1 = B(1,2)*B(0,2)*A(1,1)+B(1,2)*B(0,2)*A(1,0)+B(1,2)*B(0,2)*A(0,1)+B(1,2)*B(0,2)*A(0,0)+B(1,2)*B(0,1)*A(2,1)+B(1,2)*B(0,1)*A(2,0)+B(1,2)*B(0,1)*A(0,1)+B(1,2)*B(0,1)*A(0,0)+B(1,2)*B(0,0)*A(2,1)+B(1,2)*B(0,0)*A(2,0)+B(1,2)*B(0,0)*A(1,1)+B(1,2)*B(0,0)*A(1,0)+B(1,1)*B(0,2)*A(1,2)+B(1,1)*B(0,2)*A(1,0)+B(1,1)*B(0,2)*A(0,2)+B(1,1)*B(0,2)*A(0,0)+B(1,1)*B(0,1)*A(2,2)+B(1,1)*B(0,1)*A(2,0)+B(1,1)*B(0,1)*A(0,2)+B(1,1)*B(0,1)*A(0,0)+B(1,1)*B(0,0)*A(2,2)+B(1,1)*B(0,0)*A(2,0)+B(1,1)*B(0,0)*A(1,2)+B(1,1)*B(0,0)*A(1,0)+B(1,0)*B(0,2)*A(1,2)+B(1,0)*B(0,2)*A(1,1)+B(1,0)*B(0,2)*A(0,2)+B(1,0)*B(0,2)*A(0,1)+B(1,0)*B(0,1)*A(2,2)+B(1,0)*B(0,1)*A(2,1)+B(1,0)*B(0,1)*A(0,2)+B(1,0)*B(0,1)*A(0,1)+B(1,0)*B(0,0)*A(2,2)+B(1,0)*B(0,0)*A(2,1)+B(1,0)*B(0,0)*A(1,2)+B(1,0)*B(0,0)*A(1,1);

  A1(1,0) = A1(0,1);
  A1(2,0) = A1(0,2);
  A1(2,1) = A1(1,2);

  A2(0,0) = B(2,0)*B(0,0)*A(0,0);
  A2(0,1) = 0.5 * (B(2,1)*B(0,0)*A(0,1)+B(2,0)*B(0,1)*A(1,0));
  A2(1,1) = B(2,1)*B(0,1)*A(1,1);
  A2(0,2) = 0.5 * (B(2,2)*B(0,0)*A(0,2)+B(2,0)*B(0,2)*A(2,0));
  A2(1,2) = 0.5 * (B(2,2)*B(0,1)*A(1,2)+B(2,1)*B(0,2)*A(2,1));
  A2(2,2) = B(2,2)*B(0,2)*A(2,2);
  b2[0] = B(2,2)*B(0,0)*A(0,1)+B(2,2)*B(0,0)*A(0,0)+B(2,1)*B(0,0)*A(0,2)+B(2,1)*B(0,0)*A(0,0)+B(2,0)*B(0,2)*A(1,0)+B(2,0)*B(0,2)*A(0,0)+B(2,0)*B(0,1)*A(2,0)+B(2,0)*B(0,1)*A(0,0)+B(2,0)*B(0,0)*A(2,0)+B(2,0)*B(0,0)*A(1,0)+B(2,0)*B(0,0)*A(0,2)+B(2,0)*B(0,0)*A(0,1);
  b2[1] = B(2,2)*B(0,1)*A(1,1)+B(2,2)*B(0,1)*A(1,0)+B(2,1)*B(0,2)*A(1,1)+B(2,1)*B(0,2)*A(0,1)+B(2,1)*B(0,1)*A(2,1)+B(2,1)*B(0,1)*A(1,2)+B(2,1)*B(0,1)*A(1,0)+B(2,1)*B(0,1)*A(0,1)+B(2,1)*B(0,0)*A(2,1)+B(2,1)*B(0,0)*A(1,1)+B(2,0)*B(0,1)*A(1,2)+B(2,0)*B(0,1)*A(1,1);
  b2[2] = B(2,2)*B(0,2)*A(2,1)+B(2,2)*B(0,2)*A(2,0)+B(2,2)*B(0,2)*A(1,2)+B(2,2)*B(0,2)*A(0,2)+B(2,2)*B(0,1)*A(2,2)+B(2,2)*B(0,1)*A(0,2)+B(2,2)*B(0,0)*A(2,2)+B(2,2)*B(0,0)*A(1,2)+B(2,1)*B(0,2)*A(2,2)+B(2,1)*B(0,2)*A(2,0)+B(2,0)*B(0,2)*A(2,2)+B(2,0)*B(0,2)*A(2,1);
  c2 = B(2,2)*B(0,2)*A(1,1)+B(2,2)*B(0,2)*A(1,0)+B(2,2)*B(0,2)*A(0,1)+B(2,2)*B(0,2)*A(0,0)+B(2,2)*B(0,1)*A(2,1)+B(2,2)*B(0,1)*A(2,0)+B(2,2)*B(0,1)*A(0,1)+B(2,2)*B(0,1)*A(0,0)+B(2,2)*B(0,0)*A(2,1)+B(2,2)*B(0,0)*A(2,0)+B(2,2)*B(0,0)*A(1,1)+B(2,2)*B(0,0)*A(1,0)+B(2,1)*B(0,2)*A(1,2)+B(2,1)*B(0,2)*A(1,0)+B(2,1)*B(0,2)*A(0,2)+B(2,1)*B(0,2)*A(0,0)+B(2,1)*B(0,1)*A(2,2)+B(2,1)*B(0,1)*A(2,0)+B(2,1)*B(0,1)*A(0,2)+B(2,1)*B(0,1)*A(0,0)+B(2,1)*B(0,0)*A(2,2)+B(2,1)*B(0,0)*A(2,0)+B(2,1)*B(0,0)*A(1,2)+B(2,1)*B(0,0)*A(1,0)+B(2,0)*B(0,2)*A(1,2)+B(2,0)*B(0,2)*A(1,1)+B(2,0)*B(0,2)*A(0,2)+B(2,0)*B(0,2)*A(0,1)+B(2,0)*B(0,1)*A(2,2)+B(2,0)*B(0,1)*A(2,1)+B(2,0)*B(0,1)*A(0,2)+B(2,0)*B(0,1)*A(0,1)+B(2,0)*B(0,0)*A(2,2)+B(2,0)*B(0,0)*A(2,1)+B(2,0)*B(0,0)*A(1,2)+B(2,0)*B(0,0)*A(1,1);

  A2(1,0) = A2(0,1);
  A2(2,0) = A2(0,2);
  A2(2,1) = A2(1,2);

  A3(0,0) = B(2,0)*B(1,0)*A(0,0);
  A3(0,1) = 0.5 * (B(2,1)*B(1,0)*A(0,1)+B(2,0)*B(1,1)*A(1,0));
  A3(1,1) = B(2,1)*B(1,1)*A(1,1);
  A3(0,2) = 0.5 * (B(2,2)*B(1,0)*A(0,2)+B(2,0)*B(1,2)*A(2,0));
  A3(1,2) = 0.5 * (B(2,2)*B(1,1)*A(1,2)+B(2,1)*B(1,2)*A(2,1));
  A3(2,2) = B(2,2)*B(1,2)*A(2,2);
  b3[0] = B(2,2)*B(1,0)*A(0,1)+B(2,2)*B(1,0)*A(0,0)+B(2,1)*B(1,0)*A(0,2)+B(2,1)*B(1,0)*A(0,0)+B(2,0)*B(1,2)*A(1,0)+B(2,0)*B(1,2)*A(0,0)+B(2,0)*B(1,1)*A(2,0)+B(2,0)*B(1,1)*A(0,0)+B(2,0)*B(1,0)*A(2,0)+B(2,0)*B(1,0)*A(1,0)+B(2,0)*B(1,0)*A(0,2)+B(2,0)*B(1,0)*A(0,1);
  b3[1] = B(2,2)*B(1,1)*A(1,1)+B(2,2)*B(1,1)*A(1,0)+B(2,1)*B(1,2)*A(1,1)+B(2,1)*B(1,2)*A(0,1)+B(2,1)*B(1,1)*A(2,1)+B(2,1)*B(1,1)*A(1,2)+B(2,1)*B(1,1)*A(1,0)+B(2,1)*B(1,1)*A(0,1)+B(2,1)*B(1,0)*A(2,1)+B(2,1)*B(1,0)*A(1,1)+B(2,0)*B(1,1)*A(1,2)+B(2,0)*B(1,1)*A(1,1);
  b3[2] = B(2,2)*B(1,2)*A(2,1)+B(2,2)*B(1,2)*A(2,0)+B(2,2)*B(1,2)*A(1,2)+B(2,2)*B(1,2)*A(0,2)+B(2,2)*B(1,1)*A(2,2)+B(2,2)*B(1,1)*A(0,2)+B(2,2)*B(1,0)*A(2,2)+B(2,2)*B(1,0)*A(1,2)+B(2,1)*B(1,2)*A(2,2)+B(2,1)*B(1,2)*A(2,0)+B(2,0)*B(1,2)*A(2,2)+B(2,0)*B(1,2)*A(2,1);
  c3 = B(2,2)*B(1,2)*A(1,1)+B(2,2)*B(1,2)*A(1,0)+B(2,2)*B(1,2)*A(0,1)+B(2,2)*B(1,2)*A(0,0)+B(2,2)*B(1,1)*A(2,1)+B(2,2)*B(1,1)*A(2,0)+B(2,2)*B(1,1)*A(0,1)+B(2,2)*B(1,1)*A(0,0)+B(2,2)*B(1,0)*A(2,1)+B(2,2)*B(1,0)*A(2,0)+B(2,2)*B(1,0)*A(1,1)+B(2,2)*B(1,0)*A(1,0)+B(2,1)*B(1,2)*A(1,2)+B(2,1)*B(1,2)*A(1,0)+B(2,1)*B(1,2)*A(0,2)+B(2,1)*B(1,2)*A(0,0)+B(2,1)*B(1,1)*A(2,2)+B(2,1)*B(1,1)*A(2,0)+B(2,1)*B(1,1)*A(0,2)+B(2,1)*B(1,1)*A(0,0)+B(2,1)*B(1,0)*A(2,2)+B(2,1)*B(1,0)*A(2,0)+B(2,1)*B(1,0)*A(1,2)+B(2,1)*B(1,0)*A(1,0)+B(2,0)*B(1,2)*A(1,2)+B(2,0)*B(1,2)*A(1,1)+B(2,0)*B(1,2)*A(0,2)+B(2,0)*B(1,2)*A(0,1)+B(2,0)*B(1,1)*A(2,2)+B(2,0)*B(1,1)*A(2,1)+B(2,0)*B(1,1)*A(0,2)+B(2,0)*B(1,1)*A(0,1)+B(2,0)*B(1,0)*A(2,2)+B(2,0)*B(1,0)*A(2,1)+B(2,0)*B(1,0)*A(1,2)+B(2,0)*B(1,0)*A(1,1);

  A3(1,0) = A3(0,1);
  A3(2,0) = A3(0,2);
  A3(2,1) = A3(1,2);
}

void
polyjam::faP4P::solve( std::vector< Eigen::Vector3d> & x, std::vector< Eigen::Vector3d> & X, std::vector<Eigen::Matrix<double,3,7> > & S ) {

  Eigen::Matrix3d B = composeB(x);
  Eigen::Matrix3d G = composeG(X);
  Eigen::Matrix3d A = getA(G);

  Eigen::Matrix3d A1; Eigen::Vector3d b1; double c1;
  Eigen::Matrix3d A2; Eigen::Vector3d b2; double c2;
  Eigen::Matrix3d A3; Eigen::Vector3d b3; double c3;

  createCoeffs(A, B, A1, b1, c1, A2, b2, c2, A3, b3, c3 );
  std::vector< Eigen::Matrix<double,3,1> > solutions;
  threeConics::solve( A1, A2, A3, b1, b2, b3, c1, c2, c3, solutions );

  for( int i = 0; i < solutions.size(); i++ ) {

    Eigen::Matrix3d temp;
    temp << solutions[i][0], 1.0, 1.0,
            1.0, solutions[i][1], 1.0,
            1.0, 1.0, solutions[i][2];

    //retrieve the camera center
    Eigen::Matrix< double, 3, 1> ooo;
    ooo << 1.0,
           1.0,
           1.0;

    Eigen::Vector3d C = X[3] + G * temp.inverse() * ooo;

    //retrieve the camera matrix K
    Eigen::Matrix3d DIAC = B * temp * A * temp * B.transpose();
    Eigen::Matrix3d K = Eigen::Matrix3d::Identity();
    K(0,0) = std::sqrt(DIAC(0,0)/DIAC(2,2));
    K(1,1) = std::sqrt(DIAC(1,1)/DIAC(2,2));

    //retrieve rotation
    Eigen::Matrix3d x123; x123.col(0) = x[0]; x123.col(1) = x[1]; x123.col(2) = x[2];
    
    Eigen::Matrix3d x423; x423.col(0) = x[3]; x423.col(1) = x[1]; x423.col(2) = x[2];
    Eigen::Matrix3d x143; x143.col(0) = x[0]; x143.col(1) = x[3]; x143.col(2) = x[2];
    Eigen::Matrix3d x124; x124.col(0) = x[0]; x124.col(1) = x[1]; x124.col(2) = x[3];
    
    Eigen::Matrix3d F423; F423.col(0) = X[3]-C; F423.col(1) = X[1]-C; F423.col(2) = X[2]-C;
    Eigen::Matrix3d F143; F143.col(0) = X[0]-C; F143.col(1) = X[3]-C; F143.col(2) = X[2]-C;
    Eigen::Matrix3d F124; F124.col(0) = X[0]-C; F124.col(1) = X[1]-C; F124.col(2) = X[3]-C;

    Eigen::Matrix3d temp2;
    temp2.col(0) = (F423.determinant() / x423.determinant()) * (X[0]-C);
    temp2.col(1) = (F143.determinant() / x143.determinant()) * (X[1]-C);
    temp2.col(2) = (F124.determinant() / x124.determinant()) * (X[2]-C);

    Eigen::Matrix3d Q = K.inverse() * x123 * temp2.inverse();
    double scale = Q.col(0).norm();
    Eigen::Matrix3d R = Q / scale;
    
    Eigen::Matrix<double,3,7> Solution;
    Solution.block<3,3>(0,0) = K;
    Solution.block<3,3>(0,3) = R;
    Solution.col(6) = C;
    S.push_back(Solution);
  }
}