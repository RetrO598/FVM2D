#include <solver/FVMSolver.h>

#include <cmath>
#include <solver/numericTemplate.hpp>
#include <string>
#include <vector>

#include "pre/parameter.h"
#include "solver/turbTemplate.hpp"
#include "solver/variableDef.h"

namespace solver {
void FVMSolver::BoundaryConditions() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[ib].bnodeIndex;
    if (type == preprocess::BoundaryType::Inflow) {
      BoundInflow(ibegn, iendn);
    } else if (type == preprocess::BoundaryType::Outflow) {
      BoundOutflow(ibegn, iendn);
    } else if (type == preprocess::BoundaryType::Farfield) {
      BoundFarfield(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }

  ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    int iendn = geom.ibound[ib].bnodeIndex;
    if (param.equationtype_ == preprocess::equationType::NavierStokes &&
        (type == preprocess::BoundaryType::NoSlipWall)) {
      BoundWallVisc(ibegn, iendn);
    }

    if (param.equationtype_ == preprocess::equationType::RANS &&
        (type == preprocess::BoundaryType::NoSlipWall)) {
      BoundWallVisc(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::BoundInflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.vertexList[ib].nodeIdx;

    double ds = std::sqrt(
        geom.vertexList[ib].normal[0] * geom.vertexList[ib].normal[0] +
        geom.vertexList[ib].normal[1] * geom.vertexList[ib].normal[1]);
    double sxn = geom.vertexList[ib].normal[0] / ds;
    double syn = geom.vertexList[ib].normal[1] / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double ggm1 = dv[ibn].gamma / gam1;
    double rgas = gam1 * dv[ibn].Cp / dv[ibn].gamma;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double uabs = std::sqrt(u * u + v * v);
    double unorm = u * sxn + v * syn;
    double cosa = 0.0;
    if (uabs < 1e-16) {
      cosa = 1.0;
    } else {
      cosa = -unorm / uabs;
    }

    double c02 = dv[ibn].cs * dv[ibn].cs + 0.5 * gam1 * uabs * uabs;
    double rinv = unorm - 2.0 * dv[ibn].cs / gam1;
    double dis =
        (gam1 * cosa * cosa + 2.0) * c02 / (gam1 * rinv * rinv) - 0.5 * gam1;
    dis = std::max(dis, 1e-16);
    double cb = -rinv * (gam1 / (gam1 * cosa * cosa + 2.0)) *
                (1.0 + cosa * std::sqrt(dis));
    double cc02 = std::min((cb * cb) / c02, 1.0);
    double tb = param.TtInlet * cc02;
    double pb = param.PtInlet * std::pow((tb / param.TtInlet), ggm1);
    double rhob = pb / (rgas * tb);
    double uabsb = (2.0 * dv[ibn].Cp) * (param.TtInlet - tb);
    uabsb = std::sqrt(uabsb);
    double ub = uabsb * std::cos(param.flowAngIn);
    double vb = uabsb * std::sin(param.flowAngIn);
    CONS_VAR dummy;
    dummy.dens = rhob;
    dummy.xmom = rhob * ub;
    dummy.ymom = rhob * vb;
    dummy.ener = pb / gam1 + 0.5 * rhob * (ub * ub + vb * vb);

    auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
    auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
    CONS_VAR dummyRhs;
    numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                             dv[ibn].gamma, dv[ibn].gamma, rhs[ibn], dummyRhs);

    auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                   param.Prandtl, param.refVisc);
    if (param.equationtype_ == preprocess::equationType::NavierStokes) {
      ViscousNumericBound::ComputeResidual(
          cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
          gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
          geom.vertexList[ib].normal[1], rhs[ibn]);
    }

    if (param.equationtype_ == preprocess::equationType::RANS) {
      ViscousNumericBound::ComputeResidual(
          cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
          gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
          geom.vertexList[ib].normal[1], rhs[ibn]);

      double nuLami = dvlam[ibn].mu / cv[ibn].dens;
      double nuLamj = dummyVisc.mu / dummy.dens;

      double nuTurbi = turbVar[ibn].nu_turb;
      double nuTurbj = 3.5 * nuLamj;

      TurbSA_VAR turb_dummy;
      turb_dummy.nu_turb = nuTurbj;
      TurbSA_VAR turb_rhs_dummy;
      UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn], turb_dummy,
                                    cv[ibn], dummy, rhsTurb[ibn],
                                    turb_rhs_dummy);

      TurbViscousBound::ComputeResidual(
          nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn], gradTurbY[ibn],
          geom.vertexList[ib].normal[0], geom.vertexList[ib].normal[1],
          rhsTurb[ibn]);
    }
  }
}

void FVMSolver::BoundFarfield(int beg, int end) {
  PRIM_VAR farfield;
  if (param.farCorrect) {
    double bet = std::sqrt(1.0 - param.MaInfinity * param.MaInfinity);
    double cir = 0.25 * param.chord * cl * param.velInfinity / M_PI;
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.vertexList[ib].nodeIdx;

      double gam1 = dv[ibn].gamma - 1.0;
      double ggm1 = dv[ibn].gamma / gam1;
      double gmr = 1.0 / dv[ibn].gamma;
      double gmg = gam1 / dv[ibn].gamma;
      double xa = geom.coords[ibn].x - param.xRefPoint;
      double ya = geom.coords[ibn].y - param.yRefPoint;
      double dist = std::sqrt(xa * xa + ya * ya);
      double angle = std::atan2(ya, xa);
      double sn = std::sin(angle - param.Aoa);
      double dn = 1.0 - param.MaInfinity * param.MaInfinity * sn * sn;
      double vc = cir * bet / (dn * dist);

      farfield.velx = param.uInfinity + vc * std::sin(angle);
      farfield.vely = param.vInfinity - vc * std::cos(angle);
      double qv2 =
          farfield.velx * farfield.velx + farfield.vely * farfield.vely;
      farfield.press =
          std::pow(std::pow(param.PsInfinity, gmg) +
                       gmg * param.RhoInfinity *
                           (param.velInfinity * param.velInfinity - qv2) /
                           (2.0 * std::pow(param.PsInfinity, gmr)),
                   ggm1);
      farfield.dens =
          param.RhoInfinity * std::pow(farfield.press / param.PsInfinity, gmr);

      double ds = std::sqrt(
          geom.vertexList[ib].normal[0] * geom.vertexList[ib].normal[0] +
          geom.vertexList[ib].normal[1] * geom.vertexList[ib].normal[1]);
      double sxn = geom.vertexList[ib].normal[0] / ds;
      double syn = geom.vertexList[ib].normal[1] / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        CONS_VAR dummy;
        dummy.dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        dummy.xmom = dummy.dens * (ua + sxn * dprhoc);
        dummy.ymom = dummy.dens * (va + syn * dprhoc);
        dummy.ener = pb / gam1 +
                     0.5 * (dummy.xmom * dummy.xmom + dummy.ymom * dummy.ymom) /
                         dummy.dens;

        auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
        auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
        CONS_VAR dummyRhs;

        numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                 dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                 dummyRhs);

        auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                       param.Prandtl, param.refVisc);
        if (param.equationtype_ == preprocess::equationType::NavierStokes) {
          ViscousNumericBound::ComputeResidual(
              cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
              gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
              geom.vertexList[ib].normal[1], rhs[ibn]);
        }

        if (param.equationtype_ == preprocess::equationType::RANS) {
          ViscousNumericBound::ComputeResidual(
              cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
              gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
              geom.vertexList[ib].normal[1], rhs[ibn]);

          double nuLami = dvlam[ibn].mu / cv[ibn].dens;
          double nuLamj = dummyVisc.mu / dummy.dens;

          double nuTurbi = turbVar[ibn].nu_turb;
          double nuTurbj = 3.5 * nuLamj;

          TurbSA_VAR turb_dummy;
          turb_dummy.nu_turb = nuTurbj;
          TurbSA_VAR turb_rhs_dummy;
          UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn], turb_dummy,
                                        cv[ibn], dummy, rhsTurb[ibn],
                                        turb_rhs_dummy);

          TurbViscousBound::ComputeResidual(
              nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn], gradTurbY[ibn],
              geom.vertexList[ib].normal[0], geom.vertexList[ib].normal[1],
              rhsTurb[ibn]);
        }

      } else {
        if (qn < 0.0) {
          CONS_VAR dummy;

          dummy.dens = param.RhoInfinity;
          dummy.xmom = param.RhoInfinity * param.uInfinity;
          dummy.ymom = param.RhoInfinity * param.vInfinity;
          dummy.ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                     param.velInfinity *
                                                     param.velInfinity;
          auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
          auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
          CONS_VAR dummyRhs;
          numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                   dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                   dummyRhs);

          auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                         param.Prandtl, param.refVisc);
          if (param.equationtype_ == preprocess::equationType::NavierStokes) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);
          }

          if (param.equationtype_ == preprocess::equationType::RANS) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);

            double nuLami = dvlam[ibn].mu / cv[ibn].dens;
            double nuLamj = dummyVisc.mu / dummy.dens;

            double nuTurbi = turbVar[ibn].nu_turb;
            double nuTurbj = 3.5 * nuLamj;

            TurbSA_VAR turb_dummy;
            turb_dummy.nu_turb = nuTurbj;
            TurbSA_VAR turb_rhs_dummy;
            UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn],
                                          turb_dummy, cv[ibn], dummy,
                                          rhsTurb[ibn], turb_rhs_dummy);

            TurbViscousBound::ComputeResidual(
                nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn],
                gradTurbY[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhsTurb[ibn]);
          }
        } else {
          CONS_VAR dummy;
          dummy.dens = cv[ibn].dens;
          dummy.xmom = cv[ibn].xmom;
          dummy.ymom = cv[ibn].ymom;
          dummy.ener = cv[ibn].ener;
          auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
          auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
          CONS_VAR dummyRhs;
          numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                   dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                   dummyRhs);

          auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                         param.Prandtl, param.refVisc);
          if (param.equationtype_ == preprocess::equationType::NavierStokes) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);
          }

          if (param.equationtype_ == preprocess::equationType::RANS) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);

            double nuLami = dvlam[ibn].mu / cv[ibn].dens;
            double nuLamj = dummyVisc.mu / dummy.dens;

            double nuTurbi = turbVar[ibn].nu_turb;
            double nuTurbj = 3.5 * nuLamj;

            TurbSA_VAR turb_dummy;
            turb_dummy.nu_turb = nuTurbj;
            TurbSA_VAR turb_rhs_dummy;
            UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn],
                                          turb_dummy, cv[ibn], dummy,
                                          rhsTurb[ibn], turb_rhs_dummy);

            TurbViscousBound::ComputeResidual(
                nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn],
                gradTurbY[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhsTurb[ibn]);
          }
        }
      }
    }
  } else {
    for (int ib = beg; ib <= end; ++ib) {
      int ibn = geom.vertexList[ib].nodeIdx;

      double gam1 = dv[ibn].gamma - 1.0;
      farfield.dens = param.RhoInfinity;
      farfield.velx = param.uInfinity;
      farfield.vely = param.vInfinity;
      farfield.press = param.PsInfinity;

      double ds = std::sqrt(
          geom.vertexList[ib].normal[0] * geom.vertexList[ib].normal[0] +
          geom.vertexList[ib].normal[1] * geom.vertexList[ib].normal[1]);
      double sxn = geom.vertexList[ib].normal[0] / ds;
      double syn = geom.vertexList[ib].normal[1] / ds;

      double rhoe = cv[ibn].dens;
      double ue = cv[ibn].xmom / rhoe;
      double ve = cv[ibn].ymom / rhoe;
      double pe = dv[ibn].press;
      double qn = sxn * ue + syn * ve;
      double crho0 = dv[ibn].cs * rhoe;
      double rhoa;
      double ua;
      double va;
      double pa;
      double sgn;
      double pb;

      if (param.MaInfinity < 1.0) {
        if (qn < 0.0) {
          rhoa = farfield.dens;
          ua = farfield.velx;
          va = farfield.vely;
          pa = farfield.press;
          sgn = -1.0;
          pb = 0.5 * (pa + pe - crho0 * (sxn * (ua - ue) + syn * (va - ve)));
        } else {
          rhoa = rhoe;
          ua = ue;
          va = ve;
          pa = pe;
          sgn = 1.0;
          pb = farfield.press;
        }
        double dprhoc = sgn * (pa - pb) / crho0;
        CONS_VAR dummy;
        dummy.dens = rhoa + (pb - pa) / (dv[ibn].cs * dv[ibn].cs);
        dummy.xmom = dummy.dens * (ua + sxn * dprhoc);
        dummy.ymom = dummy.dens * (va + syn * dprhoc);
        dummy.ener = pb / gam1 +
                     0.5 * (dummy.xmom * dummy.xmom + dummy.ymom * dummy.ymom) /
                         dummy.dens;

        auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
        auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
        CONS_VAR dummyRhs;
        numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                 dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                 dummyRhs);

        auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                       param.Prandtl, param.refVisc);
        if (param.equationtype_ == preprocess::equationType::NavierStokes) {
          ViscousNumericBound::ComputeResidual(
              cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
              gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
              geom.vertexList[ib].normal[1], rhs[ibn]);
        }

        if (param.equationtype_ == preprocess::equationType::RANS) {
          ViscousNumericBound::ComputeResidual(
              cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
              gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
              geom.vertexList[ib].normal[1], rhs[ibn]);

          double nuLami = dvlam[ibn].mu / cv[ibn].dens;
          double nuLamj = dummyVisc.mu / dummy.dens;

          double nuTurbi = turbVar[ibn].nu_turb;
          double nuTurbj = 3.5 * nuLamj;

          TurbSA_VAR turb_dummy;
          turb_dummy.nu_turb = nuTurbj;
          TurbSA_VAR turb_rhs_dummy;
          UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn], turb_dummy,
                                        cv[ibn], dummy, rhsTurb[ibn],
                                        turb_rhs_dummy);

          TurbViscousBound::ComputeResidual(
              nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn], gradTurbY[ibn],
              geom.vertexList[ib].normal[0], geom.vertexList[ib].normal[1],
              rhsTurb[ibn]);
        }
      } else {
        if (qn < 0.0) {
          CONS_VAR dummy;
          dummy.dens = param.RhoInfinity;
          dummy.xmom = param.RhoInfinity * param.uInfinity;
          dummy.ymom = param.RhoInfinity * param.vInfinity;
          dummy.ener = param.PsInfinity / gam1 + 0.5 * param.RhoInfinity *
                                                     param.velInfinity *
                                                     param.velInfinity;
          auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
          auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
          CONS_VAR dummyRhs;
          numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                   dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                   dummyRhs);
          auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                         param.Prandtl, param.refVisc);
          if (param.equationtype_ == preprocess::equationType::NavierStokes) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);
          }

          if (param.equationtype_ == preprocess::equationType::RANS) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);

            double nuLami = dvlam[ibn].mu / cv[ibn].dens;
            double nuLamj = dummyVisc.mu / dummy.dens;

            double nuTurbi = turbVar[ibn].nu_turb;
            double nuTurbj = 3.5 * nuLamj;

            TurbSA_VAR turb_dummy;
            turb_dummy.nu_turb = nuTurbj;
            TurbSA_VAR turb_rhs_dummy;
            UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn],
                                          turb_dummy, cv[ibn], dummy,
                                          rhsTurb[ibn], turb_rhs_dummy);

            TurbViscousBound::ComputeResidual(
                nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn],
                gradTurbY[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhsTurb[ibn]);
          }
        } else {
          CONS_VAR dummy;
          dummy.dens = cv[ibn].dens;
          dummy.xmom = cv[ibn].xmom;
          dummy.ymom = cv[ibn].ymom;
          dummy.ener = cv[ibn].ener;
          auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
          auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
          CONS_VAR dummyRhs;
          numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                                   dv[ibn].gamma, dv[ibn].gamma, rhs[ibn],
                                   dummyRhs);

          auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                         param.Prandtl, param.refVisc);
          if (param.equationtype_ == preprocess::equationType::NavierStokes) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);
          }

          if (param.equationtype_ == preprocess::equationType::RANS) {
            ViscousNumericBound::ComputeResidual(
                cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
                gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhs[ibn]);

            double nuLami = dvlam[ibn].mu / cv[ibn].dens;
            double nuLamj = dummyVisc.mu / dummy.dens;

            double nuTurbi = turbVar[ibn].nu_turb;
            double nuTurbj = 3.5 * nuLamj;

            TurbSA_VAR turb_dummy;
            turb_dummy.nu_turb = nuTurbj;
            TurbSA_VAR turb_rhs_dummy;
            UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn],
                                          turb_dummy, cv[ibn], dummy,
                                          rhsTurb[ibn], turb_rhs_dummy);

            TurbViscousBound::ComputeResidual(
                nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn],
                gradTurbY[ibn], geom.vertexList[ib].normal[0],
                geom.vertexList[ib].normal[1], rhsTurb[ibn]);
          }
        }
      }
    }
  }
}

void FVMSolver::BoundOutflow(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    int ibn = geom.vertexList[ib].nodeIdx;

    double ds = std::sqrt(
        geom.vertexList[ib].normal[0] * geom.vertexList[ib].normal[0] +
        geom.vertexList[ib].normal[1] * geom.vertexList[ib].normal[1]);
    double sxn = geom.vertexList[ib].normal[0] / ds;
    double syn = geom.vertexList[ib].normal[1] / ds;

    double gam1 = dv[ibn].gamma - 1.0;
    double rrho = 1.0 / cv[ibn].dens;
    double u = cv[ibn].xmom * rrho;
    double v = cv[ibn].ymom * rrho;
    double q = std::sqrt(u * u + v * v);
    double mach = q / dv[ibn].cs;

    CONS_VAR dummy;
    if (mach < 1.0) {
      double rrhoc = rrho / dv[ibn].cs;
      double deltp = dv[ibn].press - param.PsOutlet;
      double rhob = cv[ibn].dens - deltp / (dv[ibn].cs * dv[ibn].cs);
      double ub = u + sxn * deltp * rrhoc;
      double vb = v + syn * deltp * rrhoc;
      double vnd = ub * sxn + vb * syn;

      if (vnd < 0.0) {
        ub = (((u) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(ub), std::abs(u));
        vb = (((v) < (0)) ? (-1.0) : (1.0)) *
             std::max(std::abs(vb), std::abs(v));
      }
      dummy.dens = rhob;
      dummy.xmom = rhob * ub;
      dummy.ymom = rhob * vb;
      dummy.ener = param.PsOutlet / gam1 + 0.5 * rhob * (ub * ub + vb * vb);
    } else {
      dummy.dens = cv[ibn].dens;
      dummy.xmom = cv[ibn].xmom;
      dummy.ymom = cv[ibn].ymom;
      dummy.ener = cv[ibn].ener;
    }

    auto dummyPrim = ConvertFromConv(dummy, dv[ibn].Cp, dv[ibn].gamma);
    auto vari = ConvertFromConv(cv[ibn], dv[ibn].Cp, dv[ibn].gamma);
    CONS_VAR dummyRhs;
    numeric->ComputeResidual(param, sxn, syn, ds, vari, dummyPrim,
                             dv[ibn].gamma, dv[ibn].gamma, rhs[ibn], dummyRhs);

    auto dummyVisc = ConvertToVISC(dummy, dv[ibn].gamma, dv[ibn].Cp,
                                   param.Prandtl, param.refVisc);
    if (param.equationtype_ == preprocess::equationType::NavierStokes) {
      ViscousNumericBound::ComputeResidual(
          cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
          gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
          geom.vertexList[ib].normal[1], rhs[ibn]);
    }

    if (param.equationtype_ == preprocess::equationType::RANS) {
      ViscousNumericBound::ComputeResidual(
          cv[ibn], dummy, dvlam[ibn], dummyVisc, gradx[ibn], grady[ibn],
          gradTx[ibn], gradTy[ibn], geom.vertexList[ib].normal[0],
          geom.vertexList[ib].normal[1], rhs[ibn]);

      double nuLami = dvlam[ibn].mu / cv[ibn].dens;
      double nuLamj = dummyVisc.mu / dummy.dens;

      double nuTurbi = turbVar[ibn].nu_turb;
      double nuTurbj = nuTurbi;

      TurbSA_VAR turb_dummy;
      turb_dummy.nu_turb = nuTurbj;
      TurbSA_VAR turb_rhs_dummy;
      UpwindTurbSA::ComputeResidual(sxn, syn, ds, turbVar[ibn], turb_dummy,
                                    cv[ibn], dummy, rhsTurb[ibn],
                                    turb_rhs_dummy);

      TurbViscousBound::ComputeResidual(
          nuLami, nuLamj, nuTurbi, nuTurbj, gradTurbX[ibn], gradTurbY[ibn],
          geom.vertexList[ib].normal[0], geom.vertexList[ib].normal[1],
          rhsTurb[ibn]);
    }
  }
}

void FVMSolver::WallVisc() {
  int ibegn = 0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (param.equationtype_ == preprocess::equationType::NavierStokes &&
        (type == preprocess::BoundaryType::NoSlipWall)) {
      BoundWallVisc(ibegn, iendn);
    }
    if (param.equationtype_ == preprocess::equationType::RANS &&
        (type == preprocess::BoundaryType::NoSlipWall)) {
      BoundWallVisc(ibegn, iendn);
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::BoundWallVisc(int beg, int end) {
  for (int ib = beg; ib <= end; ++ib) {
    // int ibn = geom.boundaryNode[ib].node;
    int ibn = geom.vertexList[ib].nodeIdx;

    cv[ibn].xmom = 0.0;
    cv[ibn].ymom = 0.0;

    if (param.equationtype_ == preprocess::equationType::RANS) {
      turbVar[ibn].nu_turb = 0.0;
    }

    ConvToDepend(ibn);
  }
}

void FVMSolver::PeriodicPrim(std::vector<PRIM_VAR> &var) {
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

        var[i].dens += var[j].dens;
        var[i].velx += var[j].velx;
        var[i].vely += var[j].vely;
        var[i].press += var[j].press;
        var[j].dens = var[i].dens;
        var[j].velx = var[i].velx;
        var[j].vely = var[i].vely;
        var[j].press = var[i].press;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicCons(std::vector<CONS_VAR> &var) {
  int ibegn = 0.0;
  for (int ib = 0; ib < geom.numBoundSegs; ++ib) {
    int iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (int ibn = ibegn; ibn <= iendn; ++ibn) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        var[i].dens += var[j].dens;
        var[i].xmom += var[j].xmom;
        var[i].ymom += var[j].ymom;
        var[i].ener += var[j].ener;
        var[j].dens = var[i].dens;
        var[j].xmom = var[i].xmom;
        var[j].ymom = var[i].ymom;
        var[j].ener = var[i].ener;
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicInt(std::vector<int> &var) {
  int i, j, ib, ibn, ibegn, iendn;

  ibegn = 0;

  for (ib = 0; ib < geom.numBoundSegs; ib++) {
    iendn = geom.ibound[ib].bnodeIndex;
    auto name = geom.bname[ib];
    auto type = param.boundaryMap.find(name)->second;
    if (type == preprocess::BoundaryType::Periodic &&
        geom.periodicMaster.find(name) != geom.periodicMaster.end()) {
      for (ibn = ibegn; ibn <= iendn; ibn++) {
        int i = geom.vertexList[ibn].nodeIdx;
        int j = geom.vertexList[ibn].periodicPair;

        var[i] += var[j];
        var[j] = var[i];
      }
    }
    ibegn = iendn + 1;
  }
}

void FVMSolver::PeriodicVisc(std::vector<double> &gradTx,
                             std::vector<double> &gradTy) {
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

        gradTx[i] += gradTx[j];
        gradTy[i] += gradTy[j];
        gradTx[j] = gradTx[i];
        gradTy[j] = gradTx[i];
      }
    }
    ibegn = iendn + 1;
  }
}
} // namespace solver