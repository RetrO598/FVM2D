#pragma once

#include <pre/geometry.h>
#include <pre/parameter.h>
#include <solver/limiter.h>
#include <solver/numeric.h>
#include <solver/variableDef.h>

#include <memory>
#include <string>
#include <vector>
#include <xtensor.hpp>
#include <xtensor/core/xtensor_forward.hpp>

// Forward-declarations for friend function signature
namespace preprocess {
template <int DIM>
struct MeshData;
}

namespace solver {
class FVMSolver;
}

namespace io {
void writeVTKFile(const std::string &filename,
                  const preprocess::MeshData<2> &mesh,
                  const solver::FVMSolver &solver);
} // namespace io

namespace solver {

class FVMSolver {
  // Grant VTK writer access to private solution data (cv, timeSteps)
  friend void io::writeVTKFile(const std::string &,
                               const preprocess::MeshData<2> &,
                               const solver::FVMSolver &);

public:
  FVMSolver(preprocess::parameter &parameter,
            const preprocess::Geometry &geometry);

  void initSolver();
  void readSolution();
  void ConvToDependAll();
  void ConvToDepend(int i);

  void BoundaryConditions();
  void BoundInflow(int beg, int end);
  void BoundOutflow(int beg, int end);
  void BoundFarfield(int beg, int end);
  void BoundWallVisc(int beg, int end);
  void WallVisc();

  void boundaryConditions();
  void boundInflow();
  void boundOutflow();
  void boundFarfield();
  void boundWallVisc();

  void Gradients();
  void GradientsVisc();

  void PeriodicCons(std::vector<CONS_VAR> &var);
  void PeriodicPrim(std::vector<PRIM_VAR> &var);
  void PeriodicInt(std::vector<int> &var);
  void PeriodicVisc(std::vector<double> &gradTx, std::vector<double> &gradTy);

  void fluxVisc();

  // void solve();

  void Timestep();

  void ZeroRes();

  void Irsmoo();

  void DensityChange(double &drho, double &drmax, int &idrmax);

  void Convergence();

  void Forces();
  void MassFlow();

  void updateResidualRK(int irk);
  void assignCVold();
  void updateCV();

  inline bool Converged() {
    if (iter >= param.maxIteration || drho <= param.convTol) {
      return true;
    }
    return false;
  }

  void writeTecplotDat();
  void writeLineDat();

  void computeWaveSpeed();
  void computeJacobianDiag(double factor);
  void lowerSweep();
  void upperSweep();
  void LUSGSupdate();
  void computeResidualLUSGS();

  // for SA Turbulent model
  void TurbGradients();
  void TurbViscous();
  void TurbSource();
  void TurbConvection();
  void TurbDissInit();
  void initTurbSolver();
  void PeriodicTurb();

  preprocess::parameter &param;
  const preprocess::Geometry &geom;
  std::unique_ptr<BaseLimiter> limiter;
  std::unique_ptr<BaseNumeric> numeric;

  int iter;

private:
  std::vector<CONS_VAR> cv;
  std::vector<CONS_VAR> cvOld;
  std::vector<CONS_VAR> diss;
  std::vector<CONS_VAR> rhs;

  std::vector<PRIM_VAR> lim;
  std::vector<PRIM_VAR> gradx;
  std::vector<PRIM_VAR> grady;
  std::vector<PRIM_VAR> umin;
  std::vector<PRIM_VAR> umax;

  std::vector<double> gradTx;
  std::vector<double> gradTy;

  std::vector<DEPEND_VAR> dv;

  std::vector<VISC_VAR> dvlam;

  // For SA turbulent model
  std::vector<TurbSA_VAR> turbVar;
  std::vector<TurbSA_VAR> turbVarOld;
  std::vector<double> gradTurbX;
  std::vector<double> gradTurbY;
  std::vector<TurbSA_VAR> dissTurb;
  std::vector<TurbSA_VAR> rhsTurb;

  std::vector<double> timeSteps;

  std::vector<CONS_VAR> rhsIter;
  std::vector<CONS_VAR> rhsOld;
  std::vector<int> nContr;

  // For LU-SGS time integrator
  std::vector<double> diag;
  std::vector<double> waveSpeed;
  std::vector<double> eigenVisc;
  std::vector<CONS_VAR> intermediateSol;
  std::vector<CONS_VAR> increment;
  std::vector<double> sourceJacobian;
  std::vector<TurbSA_VAR> intermediateSA;
  std::vector<TurbSA_VAR> incrementSA;

  double cl;
  double cd;
  double cm;
  double mflow;
  double mfratio;
  double drho;
  double drho1;
};
} // namespace solver