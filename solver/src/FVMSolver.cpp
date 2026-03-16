#include <pre/parameter.h>
#include <solver/FVMSolver.h>
#include <solver/limiter.h>
#include <solver/numeric.h>
#include <solver/variableDef.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <solver/numericFactory.hpp>
#include <vector>

#include "solver/timeIntegrator.hpp"

namespace solver {

FVMSolver::FVMSolver(preprocess::parameter &parameter,
                     const preprocess::Geometry &geometry)
    : param(parameter), geom(geometry) {
  int nNodes = geom.phyNodes;
  cv.resize(nNodes);
  cvOld.resize(nNodes);
  diss.resize(nNodes);
  rhs.resize(nNodes);

  lim.resize(nNodes);
  gradx.resize(nNodes);
  grady.resize(nNodes);
  umin.resize(nNodes);
  umax.resize(nNodes);

  timeSteps.resize(nNodes);

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    gradTx.resize(nNodes);
    gradTy.resize(nNodes);
    dvlam.resize(nNodes);
  } else if (param.equationtype_ == preprocess::equationType::RANS) {
    gradTx.resize(nNodes);
    gradTy.resize(nNodes);
    dvlam.resize(nNodes);
    turbVar.resize(nNodes);
    gradTurbX.resize(nNodes);
    gradTurbY.resize(nNodes);
    dissTurb.resize(nNodes);
    rhsTurb.resize(nNodes);
    turbVarOld.resize(nNodes);
    sourceJacobian.resize(nNodes);
    intermediateSA.resize(nNodes);
    incrementSA.resize(nNodes);
  } else {
    gradTx.resize(0);
    gradTy.resize(0);
    dvlam.resize(0);
    turbVar.resize(0);
    gradTurbX.resize(0);
    gradTurbY.resize(0);
    dissTurb.resize(0);
    rhsTurb.resize(0);
    turbVarOld.resize(0);
    sourceJacobian.resize(0);
    intermediateSA.resize(0);
    incrementSA.resize(0);
  }

  diag.resize(nNodes);
  intermediateSol.resize(nNodes);
  increment.resize(nNodes);

  dv.resize(nNodes);

  cl = 0.0;

  numeric = ConvectionFactory::create(param, geom, cv, dv, diss, rhs, lim,
                                      gradx, grady);

  limiter = LimiterFactory::create(param, geom, cv, dv, umin, umax, lim, gradx,
                                   grady);

  rhsIter.reserve(nNodes);
  rhsOld.reserve(nNodes);
  nContr.reserve(nNodes);
}

void FVMSolver::initSolver() {
  double gam1, rgas, temp, mach;

  gam1 = param.gamma - 1.0;
  rgas = gam1 * param.Cp / param.gamma;

  if (param.flowtype_ == preprocess::flowType::External) {
    CONS_VAR init;
    init.dens = param.RhoInfinity;
    init.xmom = param.RhoInfinity * param.uInfinity;
    init.ymom = param.RhoInfinity * param.vInfinity;
    init.ener = param.PsInfinity / gam1 +
                0.5 * param.RhoInfinity * param.velInfinity * param.velInfinity;
    cv.assign(geom.phyNodes, init);
  } else {
    double xmin = 1.0e+32;
    double xmax = -1.0e+32;
    for (int i = 0; i < geom.phyNodes; ++i) {
      xmin = std::min(xmin, geom.coords[i].x);
      xmax = std::max(xmax, geom.coords[i].x);
    }
    double dx = xmax - xmin;

    double pinl = param.PsRatio * param.PsOutlet;
    if (pinl >= param.PtInlet) {
      pinl = 0.99999 * param.PtInlet;
    }
    double dp = param.PsOutlet - pinl;
    double dbeta = param.approxFlowAngOut - param.flowAngIn;
    double temp =
        param.TtInlet * std::pow((pinl / param.PtInlet), (gam1 / param.gamma));
    double rho = pinl / (rgas * temp);
    double cs = std::sqrt(param.gamma * pinl / rho);
    double mach = std::sqrt(2.0 * ((param.TtInlet / temp) - 1.0) / gam1);
    double q = mach * cs;
    for (int i = 0; i < geom.phyNodes; ++i) {
      double beta = param.flowAngIn + dbeta * (geom.coords[i].x - xmin) / dx;
      double p = pinl + dp * (geom.coords[i].x - xmin) / dx;
      rho = p / (rgas * temp);
      double u = q * std::cos(beta);
      double v = q * std::sin(beta);
      cv[i].dens = rho;
      cv[i].xmom = rho * u;
      cv[i].ymom = rho * v;
      cv[i].ener = p / gam1 + 0.5 * rho * q * q;
    }
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;
        cv[i].dens = 0.5 * (cv[i].dens + cv[j].dens);
        cv[i].xmom = 0.5 * (cv[i].xmom + cv[j].xmom);
        cv[i].ymom = 0.5 * (cv[i].ymom + cv[j].ymom);
        cv[i].ener = 0.5 * (cv[i].ener + cv[j].ener);
        cv[j].dens = cv[i].dens;
        cv[j].xmom = cv[i].xmom;
        cv[j].ymom = cv[i].ymom;
        cv[j].ener = cv[i].ener;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::ConvToDependAll() {
  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;
  double g1cp = gam1 * param.Cp;
  if (param.equationtype_ == preprocess::equationType::Euler) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
      dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
      dv[i].temp = dv[i].press / (rgas * cv[i].dens);
      dv[i].cs = std::sqrt(g1cp * dv[i].temp);
      dv[i].gamma = param.gamma;
      dv[i].Cp = param.Cp;
    }
  } else {
    double s1 = 110.0;
    double s2 = 288.16;
    double s12 = 1.0 + s1 / s2;
    double cppr = param.Cp / param.Prandtl;
    for (int i = 0; i < geom.phyNodes; ++i) {
      double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
      dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
      dv[i].temp = dv[i].press / (rgas * cv[i].dens);
      dv[i].cs = std::sqrt(g1cp * dv[i].temp);
      dv[i].gamma = param.gamma;
      dv[i].Cp = param.Cp;
      double rat = std::sqrt(dv[i].temp / s2) * s12 / (1.0 + s1 / dv[i].temp);
      dvlam[i].mu = param.refVisc * rat;
      dvlam[i].lambda = dvlam[i].mu * cppr;
      if (std::isnan(rat)) {
        std::cout << geom.coords[i].x << " " << geom.coords[i].y << "\n";
        std::cout << rat << " " << dv[i].temp << " " << dv[i].press << " "
                  << param.refVisc << "\n";
        std::cout << "not a number from computing laminar mu" << "\n";
        // exit(1);
      }
      // if (param.refVisc != 0) {
      //   std::cout << param.refVisc << "\n";
      // }
    }
  }
}

void FVMSolver::ConvToDepend(int i) {
  double gam1 = param.gamma - 1.0;
  double rgas = gam1 * param.Cp / param.gamma;
  double g1cp = gam1 * param.Cp;

  if (param.equationtype_ == preprocess::equationType::Euler) {
    double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
    dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
    dv[i].temp = dv[i].press / (rgas * cv[i].dens);
    dv[i].cs = std::sqrt(g1cp * dv[i].temp);
    dv[i].gamma = param.gamma;
    dv[i].Cp = param.Cp;
  } else {
    double s1 = 110.0;
    double s2 = 288.16;
    double s12 = 1.0 + s1 / s2;

    double rhoq = cv[i].xmom * cv[i].xmom + cv[i].ymom * cv[i].ymom;
    dv[i].press = gam1 * (cv[i].ener - 0.5 * rhoq / cv[i].dens);
    dv[i].temp = dv[i].press / (rgas * cv[i].dens);
    dv[i].cs = std::sqrt(g1cp * dv[i].temp);
    dv[i].gamma = param.gamma;
    dv[i].Cp = param.Cp;
    double rat = std::sqrt(dv[i].temp / s2) * s12 / (1.0 + s1 / dv[i].temp);
    dvlam[i].mu = param.refVisc * rat;
    dvlam[i].lambda = dvlam[i].mu * (param.Cp / param.Prandtl);
  }
}

void FVMSolver::updateResidualRK(int irk) {
  numeric->DissipInit();
  if (param.equationtype_ == preprocess::equationType::RANS) {
    TurbDissInit();
  }

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    GradientsVisc();
    fluxVisc();
  } else if (param.equationtype_ == preprocess::equationType::Euler) {
    Gradients();
  } else if (param.equationtype_ == preprocess::equationType::RANS) {
    GradientsVisc();
    fluxVisc();
    TurbGradients();
    TurbViscous();
  }

  limiter->limiterInit();
  limiter->limiterUpdate();

  numeric->FluxNumeric();
  if (param.equationtype_ == preprocess::equationType::RANS) {
    TurbConvection();
    TurbSource();
  }

  BoundaryConditions();
  ZeroRes();
  PeriodicCons(rhs);
  if (param.equationtype_ == preprocess::equationType::RANS) {
    PeriodicTurb();
  }

  double fac = param.stageCoeff[irk] * param.CFL;
  for (int i = 0; i < geom.phyNodes; ++i) {
    double adtv = fac * timeSteps[i] / geom.vol[i];
    rhs[i].dens *= adtv;
    rhs[i].xmom *= adtv;
    rhs[i].ymom *= adtv;
    rhs[i].ener *= adtv;
  }

  if (param.equationtype_ == preprocess::equationType::RANS) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      double adtv = fac * timeSteps[i] / geom.vol[i];
      rhsTurb[i].nu_turb *= adtv;
    }
  }
}

void FVMSolver::assignCVold() {
  cvOld = cv;
  if (param.equationtype_ == preprocess::equationType::RANS) {
    turbVarOld = turbVar;
  }
}

void FVMSolver::updateCV() {
  for (int i = 0; i < geom.phyNodes; ++i) {
    cv[i].dens = cvOld[i].dens - rhs[i].dens;
    cv[i].xmom = cvOld[i].xmom - rhs[i].xmom;
    cv[i].ymom = cvOld[i].ymom - rhs[i].ymom;
    cv[i].ener = cvOld[i].ener - rhs[i].ener;
    if (std::isnan(rhs[i].ener)) {
      std::cout << "not a number from update energy" << "\n";
    }
  }

  if (param.equationtype_ == preprocess::equationType::RANS) {
    for (int i = 0; i < geom.phyNodes; ++i) {
      turbVar[i].nu_turb =
          std::max(turbVarOld[i].nu_turb - rhsTurb[i].nu_turb, 0.0);
      if (std::isnan(turbVar[i].nu_turb)) {
        std::cout << "not a number from update" << "\n";
      }
    }
  }
}

void FVMSolver::computeWaveSpeed() {
  std::size_t sum = 0;
  for (auto &i : geom.pointList) {
    sum += i.getnPoint();
  }

  waveSpeed.resize(sum);
  eigenVisc.resize(sum);

  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    for (std::size_t j = 0; j < geom.pointList[iPoint].getnPoint(); ++j) {
      auto jPoint = geom.pointList[iPoint].getPoint(j);
      auto jEdge = geom.pointList[iPoint].getEdge(j);

      double sx = geom.sij[jEdge].x;
      double sy = geom.sij[jEdge].y;

      double ds = std::sqrt(sx * sx + sy * sy);

      double nx = sx / ds;
      double ny = sy / ds;

      double velx = 0.5 * (cv[iPoint].xmom / cv[iPoint].dens +
                           cv[jPoint].xmom / cv[jPoint].dens);
      double vely = 0.5 * (cv[iPoint].ymom / cv[iPoint].dens +
                           cv[jPoint].ymom / cv[jPoint].dens);

      double cs = 0.5 * (dv[iPoint].cs + dv[jPoint].cs);

      double wavespeed = (std::abs(velx * nx + vely * ny) + cs) * ds;

      double uj = cv[jPoint].xmom / cv[jPoint].dens;
      double vj = cv[jPoint].ymom / cv[jPoint].dens;

      waveSpeed[index] = wavespeed;

      if (param.equationtype_ == preprocess::equationType::NavierStokes) {
        double rhoAve = 0.5 * (cv[iPoint].dens + cv[jPoint].dens);
        double gammaAve = 0.5 * (dv[iPoint].gamma + dv[jPoint].gamma);
        double mue = 0.5 * (dvlam[iPoint].mu / param.Prandtl +
                            dvlam[jPoint].mu / param.Prandtl);
        double f1 = std::max(4.0 / 3.0 / rhoAve, gammaAve / rhoAve);
        double eigen = f1 * mue;
        eigenVisc[index] = eigen * ds * ds;
      }

      if (param.equationtype_ == preprocess::equationType::RANS) {
        double muTurbi = 0.0;
        double muTurbj = 0.0;
        double lambdaTurbi = 0.0;
        double lambdaTurbj = 0.0;

        double nuLami = dvlam[iPoint].mu / cv[iPoint].dens;
        double nuLamj = dvlam[jPoint].mu / cv[jPoint].dens;
        double cv1_3 = 357.911;

        double Ji_i = turbVar[iPoint].nu_turb / nuLami;
        double Ji2_i = Ji_i * Ji_i;
        double Ji3_i = Ji2_i * Ji_i;
        double fv1_i = Ji3_i / (Ji3_i + cv1_3);

        double Ji_j = turbVar[jPoint].nu_turb / nuLamj;
        double Ji2_j = Ji_j * Ji_j;
        double Ji3_j = Ji2_j * Ji_j;
        double fv1_j = Ji3_j / (Ji3_j + cv1_3);
        muTurbi = turbVar[iPoint].nu_turb * cv[iPoint].dens * fv1_i;
        muTurbj = turbVar[jPoint].nu_turb * cv[jPoint].dens * fv1_j;

        double rhoAve = 0.5 * (cv[iPoint].dens + cv[jPoint].dens);
        double gammaAve = 0.5 * (dv[iPoint].gamma + dv[jPoint].gamma);
        double mue = 0.5 * (dvlam[iPoint].mu / param.Prandtl +
                            dvlam[jPoint].mu / param.Prandtl + muTurbi / 0.9 +
                            muTurbj / 0.9);
        double f1 = std::max(4.0 / 3.0 / rhoAve, gammaAve / rhoAve);
        double eigen = f1 * mue;
        eigenVisc[index] = eigen * ds * ds;
      }
      index++;
    }
  }
}

void FVMSolver::computeJacobianDiag(double factor) {
  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    diag[iPoint] = geom.vol[iPoint] / (timeSteps[iPoint] * param.CFL);
    for (std::size_t j = 0; j < geom.pointList[iPoint].getnPoint(); ++j) {
      double lambda = waveSpeed[index] * factor;
      waveSpeed[index] = lambda;
      diag[iPoint] += 0.5 * lambda;
      if (param.equationtype_ == preprocess::equationType::NavierStokes ||
          param.equationtype_ == preprocess::equationType::RANS) {
        diag[iPoint] += eigenVisc[index] / geom.vol[iPoint];
      }
      index++;
    }
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[ib].bnodeIndex;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        diag[i] += diag[j];
        diag[j] = diag[i];
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::lowerSweep() {
  CONS_VAR df{0.0, 0.0, 0.0, 0.0};
  TurbSA_VAR dfSA{0.0};
  std::size_t index = 0;
  for (std::size_t iPoint = 0; iPoint < geom.pointList.size(); ++iPoint) {
    df = {0.0, 0.0, 0.0, 0.0};
    dfSA = {0.0};
    for (std::size_t ineighbor = 0;
         ineighbor < geom.pointList[iPoint].getnPoint(); ++ineighbor) {
      std::size_t jPoint = geom.pointList[iPoint].getPoint(ineighbor);
      std::size_t jEdge = geom.pointList[iPoint].getEdge(ineighbor);

      if (jPoint < iPoint) {
        CONS_VAR sol = cv[jPoint];
        CONS_VAR dSol;
        dSol.dens = sol.dens + intermediateSol[jPoint].dens;
        dSol.xmom = sol.xmom + intermediateSol[jPoint].xmom;
        dSol.ymom = sol.ymom + intermediateSol[jPoint].ymom;
        dSol.ener = sol.ener + intermediateSol[jPoint].ener;

        double sx, sy, nx, ny, gamma;

        sx = geom.sij[jEdge].x;
        sy = geom.sij[jEdge].y;
        double ds = std::sqrt(sx * sx + sy * sy);
        nx = sx / ds;
        ny = sy / ds;

        if (geom.edge[jEdge].nodei == jPoint) {
          nx = -nx;
          ny = -ny;
        }

        gamma = dv[jPoint].gamma;

        TurbSA_VAR solTurb = turbVar[jPoint];
        TurbSA_VAR dSolTurb;
        dSolTurb.nu_turb =
            turbVar[jPoint].nu_turb + intermediateSA[jPoint].nu_turb;

        TurbSA_VAR fTurb;
        TurbSA_VAR dfTurb;
        if (param.equationtype_ == preprocess::equationType::RANS) {
          fTurb = turbFluxDifference::computeDifference(sol, solTurb, nx, ny);
          dfTurb =
              turbFluxDifference::computeDifference(dSol, dSolTurb, nx, ny);

          dfTurb.nu_turb -= fTurb.nu_turb;
        }

        auto f = fluxDifference::computeDifference(sol, nx, ny, gamma);
        auto dfj = fluxDifference::computeDifference(dSol, nx, ny, gamma);

        dfj.dens -= f.dens;
        dfj.xmom -= f.xmom;
        dfj.ymom -= f.ymom;
        dfj.ener -= f.ener;

        double ra = waveSpeed[index];
        if (param.equationtype_ == preprocess::equationType::NavierStokes ||
            param.equationtype_ == preprocess::equationType::RANS) {
          double iPointx = geom.pointList[iPoint].getCoord(0);
          double iPointy = geom.pointList[iPoint].getCoord(1);
          double jPointx = geom.pointList[jPoint].getCoord(0);
          double jPointy = geom.pointList[jPoint].getCoord(1);
          double length = std::sqrt((jPointx - iPointx) * (jPointx - iPointx) +
                                    (jPointy - iPointy) * (jPointy - iPointy));
          ra += eigenVisc[index] / ds / length;
        }

        df.dens += (dfj.dens * ds - ra * intermediateSol[jPoint].dens);
        df.xmom += (dfj.xmom * ds - ra * intermediateSol[jPoint].xmom);
        df.ymom += (dfj.ymom * ds - ra * intermediateSol[jPoint].ymom);
        df.ener += (dfj.ener * ds - ra * intermediateSol[jPoint].ener);

        if (param.equationtype_ == preprocess::equationType::RANS) {
          dfSA.nu_turb +=
              (dfTurb.nu_turb * ds - ra * intermediateSA[jPoint].nu_turb);
        }
      }
      index++;
    }

    intermediateSol[iPoint].dens =
        (-rhs[iPoint].dens - 0.5 * df.dens) / diag[iPoint];
    intermediateSol[iPoint].xmom =
        (-rhs[iPoint].xmom - 0.5 * df.xmom) / diag[iPoint];
    intermediateSol[iPoint].ymom =
        (-rhs[iPoint].ymom - 0.5 * df.ymom) / diag[iPoint];
    intermediateSol[iPoint].ener =
        (-rhs[iPoint].ener - 0.5 * df.ener) / diag[iPoint];

    if (param.equationtype_ == preprocess::equationType::RANS) {
      intermediateSA[iPoint].nu_turb =
          (-rhsTurb[iPoint].nu_turb - 0.5 * dfSA.nu_turb) /
          (diag[iPoint] - sourceJacobian[iPoint]);
    }
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[ib].bnodeIndex;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        intermediateSol[i].dens +=
            (intermediateSol[j].dens + rhs[j].dens / diag[j]);
        intermediateSol[i].xmom +=
            (intermediateSol[j].xmom + rhs[j].xmom / diag[j]);
        intermediateSol[i].ymom +=
            (intermediateSol[j].ymom + rhs[j].ymom / diag[j]);
        intermediateSol[i].ener +=
            (intermediateSol[j].ener + rhs[j].ener / diag[j]);
        intermediateSol[j].dens = intermediateSol[i].dens;
        intermediateSol[j].xmom = intermediateSol[i].xmom;
        intermediateSol[j].ymom = intermediateSol[i].ymom;
        intermediateSol[j].ener = intermediateSol[i].ener;

        if (param.equationtype_ == preprocess::equationType::RANS) {
          intermediateSA[i].nu_turb +=
              (intermediateSA[j].nu_turb +
               rhsTurb[j].nu_turb / (diag[j] - sourceJacobian[j]));
          intermediateSA[j].nu_turb = intermediateSA[i].nu_turb;
        }
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::upperSweep() {
  CONS_VAR df{0.0, 0.0, 0.0, 0.0};
  TurbSA_VAR dfSA{0.0};
  std::size_t index = waveSpeed.size();
  for (int iPoint = geom.pointList.size() - 1; iPoint > -1; --iPoint) {
    df = {0.0, 0.0, 0.0, 0.0};
    dfSA = {0.0};
    for (int ineighbor = geom.pointList[iPoint].getnPoint() - 1; ineighbor > -1;
         --ineighbor) {
      index--;
      std::size_t jPoint = geom.pointList[iPoint].getPoint(ineighbor);
      std::size_t jEdge = geom.pointList[iPoint].getEdge(ineighbor);

      if (jPoint > iPoint) {
        CONS_VAR sol = cv[jPoint];
        CONS_VAR dSol;
        dSol.dens = sol.dens + increment[jPoint].dens;
        dSol.xmom = sol.xmom + increment[jPoint].xmom;
        dSol.ymom = sol.ymom + increment[jPoint].ymom;
        dSol.ener = sol.ener + increment[jPoint].ener;

        double sx, sy, nx, ny, gamma;

        sx = geom.sij[jEdge].x;
        sy = geom.sij[jEdge].y;
        double ds = std::sqrt(sx * sx + sy * sy);
        nx = sx / ds;
        ny = sy / ds;

        if (geom.edge[jEdge].nodei == jPoint) {
          nx = -nx;
          ny = -ny;
        }

        gamma = dv[jPoint].gamma;

        TurbSA_VAR solTurb = turbVar[jPoint];
        TurbSA_VAR dSolTurb;
        dSolTurb.nu_turb =
            turbVar[jPoint].nu_turb + incrementSA[jPoint].nu_turb;

        TurbSA_VAR fTurb;
        TurbSA_VAR dfTurb;
        if (param.equationtype_ == preprocess::equationType::RANS) {
          fTurb = turbFluxDifference::computeDifference(sol, solTurb, nx, ny);
          dfTurb =
              turbFluxDifference::computeDifference(dSol, dSolTurb, nx, ny);

          dfTurb.nu_turb -= fTurb.nu_turb;
        }

        auto f = fluxDifference::computeDifference(sol, nx, ny, gamma);
        auto dfj = fluxDifference::computeDifference(dSol, nx, ny, gamma);

        dfj.dens -= f.dens;
        dfj.xmom -= f.xmom;
        dfj.ymom -= f.ymom;
        dfj.ener -= f.ener;

        double ra = waveSpeed[index];
        if (param.equationtype_ == preprocess::equationType::NavierStokes ||
            param.equationtype_ == preprocess::equationType::RANS) {
          double iPointx = geom.pointList[iPoint].getCoord(0);
          double iPointy = geom.pointList[iPoint].getCoord(1);
          double jPointx = geom.pointList[jPoint].getCoord(0);
          double jPointy = geom.pointList[jPoint].getCoord(1);
          double length = std::sqrt((jPointx - iPointx) * (jPointx - iPointx) +
                                    (jPointy - iPointy) * (jPointy - iPointy));
          ra += eigenVisc[index] / ds / length;
        }

        df.dens += (dfj.dens * ds - ra * increment[jPoint].dens);
        df.xmom += (dfj.xmom * ds - ra * increment[jPoint].xmom);
        df.ymom += (dfj.ymom * ds - ra * increment[jPoint].ymom);
        df.ener += (dfj.ener * ds - ra * increment[jPoint].ener);

        if (param.equationtype_ == preprocess::equationType::RANS) {
          dfSA.nu_turb +=
              (dfTurb.nu_turb * ds - ra * incrementSA[jPoint].nu_turb);
        }
      }
    }

    increment[iPoint].dens =
        intermediateSol[iPoint].dens - 0.5 * df.dens / diag[iPoint];
    increment[iPoint].xmom =
        intermediateSol[iPoint].xmom - 0.5 * df.xmom / diag[iPoint];
    increment[iPoint].ymom =
        intermediateSol[iPoint].ymom - 0.5 * df.ymom / diag[iPoint];
    increment[iPoint].ener =
        intermediateSol[iPoint].ener - 0.5 * df.ener / diag[iPoint];

    if (param.equationtype_ == preprocess::equationType::RANS) {
      incrementSA[iPoint].nu_turb =
          intermediateSA[iPoint].nu_turb -
          0.5 * dfSA.nu_turb / (diag[iPoint] - sourceJacobian[iPoint]);
    }
  }

  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[ib].bnodeIndex;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        increment[i].dens +=
            (increment[j].dens - intermediateSol[j].dens / diag[j]);
        increment[i].xmom +=
            (increment[j].xmom - intermediateSol[j].xmom / diag[j]);
        increment[i].ymom +=
            (increment[j].ymom - intermediateSol[j].ymom / diag[j]);
        increment[i].ener +=
            (increment[j].ener - intermediateSol[j].ener / diag[j]);
        increment[j].dens = increment[i].dens;
        increment[j].xmom = increment[i].xmom;
        increment[j].ymom = increment[i].ymom;
        increment[j].ener = increment[i].ener;

        if (param.equationtype_ == preprocess::equationType::RANS) {
          incrementSA[i].nu_turb +=
              (incrementSA[j].nu_turb -
               intermediateSA[j].nu_turb / (diag[j] - sourceJacobian[j]));
          incrementSA[j].nu_turb = incrementSA[i].nu_turb;
        }
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::LUSGSupdate() {
  for (std::size_t i = 0; i < geom.phyNodes; ++i) {
    cv[i].dens += increment[i].dens;
    cv[i].xmom += increment[i].xmom;
    cv[i].ymom += increment[i].ymom;
    cv[i].ener += increment[i].ener;
  }

  if (param.equationtype_ == preprocess::equationType::RANS) {
    for (std::size_t i = 0; i < geom.phyNodes; ++i) {
      turbVar[i].nu_turb += std::max(incrementSA[i].nu_turb, 0.0);
    }
  }
}

void FVMSolver::computeResidualLUSGS() {
  numeric->DissipInit();
  if (param.equationtype_ == preprocess::equationType::RANS) {
    TurbDissInit();
  }

  if (param.equationtype_ == preprocess::equationType::NavierStokes) {
    GradientsVisc();
    fluxVisc();
  } else if (param.equationtype_ == preprocess::equationType::Euler) {
    Gradients();
  } else if (param.equationtype_ == preprocess::equationType::RANS) {
    GradientsVisc();
    fluxVisc();
    TurbGradients();
    TurbViscous();
  }

  limiter->limiterInit();
  limiter->limiterUpdate();

  numeric->FluxNumeric();
  if (param.equationtype_ == preprocess::equationType::RANS) {
    TurbConvection();
    TurbSource();
  }

  BoundaryConditions();
  ZeroRes();
  PeriodicCons(rhs);
  if (param.equationtype_ == preprocess::equationType::RANS) {
    PeriodicTurb();
  }
}
} // namespace solver