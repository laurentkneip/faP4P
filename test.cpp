#include <iostream>
#include <sys/time.h>

#include "faP4P.hpp"

using namespace std;
using namespace Eigen;

Vector3d
generateRandomPoint( double maximumDepth, double minimumDepth )
{
  Vector3d cleanPoint;
  cleanPoint[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  cleanPoint[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  Vector3d direction = cleanPoint / cleanPoint.norm();
  cleanPoint =
      (maximumDepth-minimumDepth) * cleanPoint + minimumDepth * direction;
  return cleanPoint;
}

Vector3d
generateRandomTranslation( double maximumParallax )
{
  Vector3d translation;
  translation[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  translation[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
  return maximumParallax * translation;
}

Matrix3d
generateRandomRotation( double maxAngle )
{
  Vector3d rpy;
  rpy[0] = ((double) rand())/ ((double) RAND_MAX);
  rpy[1] = ((double) rand())/ ((double) RAND_MAX);
  rpy[2] = ((double) rand())/ ((double) RAND_MAX);

  rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
  rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
  rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

  Matrix3d R1;
  R1(0,0) = 1.0;
  R1(0,1) = 0.0;
  R1(0,2) = 0.0;
  R1(1,0) = 0.0;
  R1(1,1) = cos(rpy[0]);
  R1(1,2) = -sin(rpy[0]);
  R1(2,0) = 0.0;
  R1(2,1) = -R1(1,2);
  R1(2,2) = R1(1,1);

  Matrix3d R2;
  R2(0,0) = cos(rpy[1]);
  R2(0,1) = 0.0;
  R2(0,2) = sin(rpy[1]);
  R2(1,0) = 0.0;
  R2(1,1) = 1.0;
  R2(1,2) = 0.0;
  R2(2,0) = -R2(0,2);
  R2(2,1) = 0.0;
  R2(2,2) = R2(0,0);

  Matrix3d R3;
  R3(0,0) = cos(rpy[2]);
  R3(0,1) = -sin(rpy[2]);
  R3(0,2) = 0.0;
  R3(1,0) =-R3(0,1);
  R3(1,1) = R3(0,0);
  R3(1,2) = 0.0;
  R3(2,0) = 0.0;
  R3(2,1) = 0.0;
  R3(2,2) = 1.0;

  Matrix3d rotation = R3 * R2 * R1;

  rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
  rotation.col(2) = rotation.col(0).cross(rotation.col(1));
  rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
  rotation.col(1) = rotation.col(2).cross(rotation.col(0));
  rotation.col(1) = rotation.col(1) / rotation.col(1).norm();

  return rotation;
}


int main( int argc, char** argv )
{
  //initialize random seed
  struct timeval tic;
  gettimeofday( &tic, 0 );
  srand ( tic.tv_usec );
  
  //generate random points and then move forward along principal axis so they lie in front of the camera
  vector< Vector3d > X;
  for( int i = 0; i < 4; i++ ) {

    X.push_back(generateRandomPoint(0.0,4.0));
    X.back()[2] += 8.0;

  }

  //invent small translations and rotations (to make sure the points stay in the field of view)
  Vector3d C = generateRandomTranslation( 0.1 );
  Matrix3d R = generateRandomRotation( 0.1 );

  //invent some camera matrix
  Matrix3d K;
  K << 400.0,   0.0, 350.0,
         0.0, 300.0, 350.0,
         0.0,   0.0,   1.0;

  //calculate the reprojections into the image plane
  vector< Vector3d > f;
  for ( int i = 0; i < 4; i++ ) {

    Vector3d temp = K * R * (X[i] - C);
    double scale = temp[2];
    Vector3d temp2 = temp / scale;
    f.push_back(temp2);
  }

  //calculate the image points normalized by the principal point and some approximate focal length (serves to normalize algorithm input)
  Matrix3d Kapprox;
  Kapprox << 350.0,   0.0, 350.0,
               0.0, 350.0, 350.0,
               0.0,   0.0,   1.0;
  vector< Vector3d > x;
  for ( int i = 0; i < 4; i++ )
    x.push_back( Kapprox.inverse() * f[i] );

  //print the groundtruth variables
  cout << "Groundtruth camera matrix:\n" << (Kapprox.inverse() * K) << "\n";
  cout << "Groundtruth camera orientation:\n" << R << "\n";
  cout << "Groundtruth camera center:\n"  << C.transpose() << "\n\n";
  
  //call the solver (results are now internally printed, feel free to change that)
  vector<Matrix<double,3,7> > S;
  polyjam::faP4P::solve(x,X,S);

  //print the result
  cout << "Solver found " << S.size() << " solutions. Printing them:\n\n";
  for( int i = 0; i < S.size(); i++ ) {
    cout << "Solution " << i << ":\n";
    cout << "Camera matrix K:\n" << S[i].block<3,3>(0,0) << "\n";
    cout << "Camera rotation R:\n" << S[i].block<3,3>(0,3) << "\n";
    cout << "Camera center C:\n" << S[i].col(6).transpose() << "\n\n";
  }

  return 0;
}