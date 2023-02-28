#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include "./Eigen/Sparse" 


// Element force vector
void ComputeElementForceVector(const double& Xa, const double& Xb, double* frc_vec);

// Element mass matrix
void ComputeElementMassMatrix(const double& Xa, const double& Xb, double mass_mat[2][2]);

// Compute element's contribution to the L2 error = Integral from Xa to Xb of (f-Phf)^2.
double ComputeElementL2Error(const double& Xa, const double& Xb, const double& Phf_a, const double& Phf_b);

int main()
{
  // Create a 1D grid over the interval [0,1]
  // Set the grid size over [0,0.5] = h
  // Set the grid size over [0.5,1] = 2h
  const int nElements = 1;
  const int nNodes = nElements+1;
  std::vector<double> coordinates(nNodes);
  for(int n=0; n<nNodes; ++n)
    {
      // Q1. INITIALIZE NODAL COORDINATES
    }
  
  // Global load vector
  Eigen::VectorXd loadVec(nNodes);
  // Q2. INITIALIZE LOAD VECTOR TO ZERO 

  // Entries of the global mass matrix
  std::vector<Eigen::Triplet<double>> matTriplets;
  matTriplets.reserve(nElements*4);

  // Compute entries of the mass matrix and load vector
  double elm_frc_vec[2];
  double elm_mass_mat[2][2];
  for(int e=0; e<nElements; ++e)
    {
      const int a = e;
      const int b = e+1;

      const double& Xa = coordinates[a];
      const double& Xb = coordinates[b];

      // compute the element force vector and assemble it into the global vector
      ComputeElementForceVector(Xa, Xb, elm_frc_vec);

      // Q3. ASSEMBLE ELEMENT FORCE VECTOR ENTRIES INTO "loadVec"
      // ...
      
      // compute the element mass matrix
      ComputeElementMassMatrix(Xa, Xb, elm_mass_mat);

      // Q4. APPEND ELEMENT MASS MATRIX ENTRIES TO "matTriplets"
      // ....
    }
  matTriplets.shrink_to_fit();
  
  // Create the global mass matrix and assemble it
  Eigen::SparseMatrix<double, Eigen::RowMajor> massMat(nNodes, nNodes);
  massMat.setFromTriplets(matTriplets.begin(), matTriplets.end());

  // Linear solver
  Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
  solver.analyzePattern(massMat);
  solver.factorize(massMat);
  
  // Dof vector
  Eigen::VectorXd dofVec(nNodes);

  // Solve
  dofVec = solver.solve(loadVec);

  // Check that the solve went well
  assert(solver.info()==Eigen::Success && "Linear solve was unsuccessful");

  // Visualize the solution
  std::fstream file;
  file.open("sol.dat", std::ios::out);
  file<<"x \t Phf(x)";
  for(int n=0; n<nNodes; ++n)
    file <<"\n"<<coordinates[n]<<"\t"<<dofVec(n);
  file.close();

  // Compute the L2 norm of the error
  double L2_error = 0.;
  for(int e=0; e<nElements; ++e)
    {
      const double& Xa = coordinates[e];
      const double& Xb = coordinates[e+1];
      const double& Phf_a = dofVec(e);
      const double& Phf_b = dofVec(e+1);
      
      L2_error += ComputeElementL2Error(Xa, Xb, Phf_a, Phf_b);
    }
  L2_error = std::sqrt(L2_error);
  std::cout << std::endl << "L2 error: " << L2_error << std::endl;
  
  // done
}


// Element force vector
void ComputeElementForceVector(const double& Xa, const double& Xb, double* frc_vec)
{
  // Q5 EVALUATE THE ELEMENT FORCE VECTOR TO "frc_vec"
  // EVALUATE INTEGRALS EXACTLY
  // ..

  // done
  return;
}


// Element mass matrix
void ComputeElementMassMatrix(const double& Xa, const double& Xb,
			      double mass_mat[2][2])
{
  // Q6 EVALUATE THE ELEMENT MASS MATRIX TO "mass_mat"
  // EVALUATE INTEGRALS EXACTLY
  // ..

  // done
  return;
}


// // Compute element's contribution to the L2 error = Integral from Xa to Xb of (f-Phf)^2.
double ComputeElementL2Error(const double& Xa, const double& Xb, const double& Phf_a, const double& Phf_b)
{
  double integral = 0.;
  
  // Q7 INTEGRATE (f-Phf)^2 FROM Xa TO Xb
  // EVALUATE INTEGRALS EXACTLY
  // ..

  // done
  return integral;
}

