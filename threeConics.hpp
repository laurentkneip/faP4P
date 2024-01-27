/*************************************************************************
 *                                                                       *
 * polyjam, a polynomial solver generator for C++                        *
 * Copyright (C) 2015 Laurent Kneip, The Australian National University  *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *                                                                       *
 *************************************************************************/

//This code is automatically generated by polyjam for solving threeConics.
//It is licensed under the GNU GPL terms.
//Please contact the author of polyjam for proprietary use.

#ifndef POLYJAM_threeConics_HPP_
#define POLYJAM_threeConics_HPP_

#include <stdlib.h>
#include <Eigen/Eigen>
#include <vector>
#include <list>


namespace polyjam
{
namespace threeConics
{

  void initRow(
      Eigen::MatrixXd  & M2,
      const Eigen::MatrixXd & M1,
      int row2,
      int row1,
      const int * cols2,
      const int * cols1,
      size_t numberCols );

  void solve( Eigen::Matrix3d & A1, Eigen::Matrix3d & A2, Eigen::Matrix3d & A3, Eigen::Vector3d & b1, Eigen::Vector3d & b2, Eigen::Vector3d & b3, double c1, double c2, double c3, std::vector< Eigen::Matrix<double,3,1> > & solutions );

}
}

#endif /* POLYJAM_threeConics_HPP_ */