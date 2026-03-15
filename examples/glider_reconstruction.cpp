#include "matplot/freestanding/axes_functions.h"
#include "unscented/primitives.h"
#include "unscented/ukf.h"

#include <matplot/matplot.h>
#include <iostream>
#include <random>

using Pos = unscented::Scalar;
using Vel = unscented::Scalar;

//////////////////////////////////////////////////////////////////////////////
// Set up the state
//////////////////////////////////////////////////////////////////////////////

enum StateElements
{
  POS_X = 0,
  POS_Y,
  POS_Z,
  VEL_U,
  VEL_W,
  AVEL_Q,
  QUATERNION,
};

using State =
    unscented::Compound<unscented::Vector<6>, unscented::UnitQuaternion>;

using DerivState = unscented::Vector<10>;

using Vec3 = std::tuple<double, double, double>;

// x,y, z, u, w, q, q0, q1, q2, q3

//////////////////////////////////////////////////////////////////////////////
// Set up the system model
//////////////////////////////////////////////////////////////////////////////

struct GliderAero
{
  double S = 10.0;
  double dcldalpha = 6.0;
  double k = 0.00769772063821893;
  double g = 9.81;
  double CL0 = 0.5248831739026187;
  double CD0 = 0.0107362467389172;
  double m = 500.0;
  double rho = 1.225; // air density kg/m^3
};

/*
{'LD_best': 55,
 'V_LDbest': 30.8641975308642,
 'V_cruise': 20.5761316872428,
 'S': 10.0, 'm': 500.0, 'rho': 1.225, 'g': 9.881,
 'dcldalpha': 6.0,
 'wind_n': 2.5308641975308643, 'wind_e': -11.882716049382717,
 'k': 0.00769772063821893, 'CD0': 0.0107362467389172, 'CL0': 0.5248831739026187}
*/

GliderAero parms;

DerivState system_ode(const State& state)
{
  auto& [states, attitude] = state.data;
  auto& u = states[VEL_U];
  auto& w = states[VEL_W];
  auto& q = states[AVEL_Q];
  auto& quaternion = attitude.get_q().coeffs();
  auto& q0 = quaternion[0];
  auto& q1 = quaternion[1];
  auto& q2 = quaternion[2];
  auto& q3 = quaternion[3];

  auto V_sq = u * u + w * w;
  auto V = sqrt(V_sq);
  auto alpha = atan(w / u);

  auto CL = parms.dcldalpha * alpha + parms.CL0;
  // CL = sympy.Min(1.5, sympy.Max(0.0, cl_raw))
  auto CD = parms.CD0 + parms.k * CL * CL;
  // CLL = sympy.Min(CL, 1.0)
  auto CLL = CL;

  auto Q = 0.5 * parms.rho * V_sq;
  auto QS_m = Q * parms.S / parms.m;
  auto sa = w / V; // sin(alpha)
  auto ca = u / V; // cos(alpha)
  auto Cx = CLL * sa - CD * ca;
  auto Cz = -CLL * ca - CD * sa;

  // orientation
  auto a1 = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
  auto b1 = 2 * (q1 * q2 - q0 * q3);
  auto c1 = 2 * (q0 * q2 + q1 * q3);

  auto a2 = 2 * (q1 * q2 + q0 * q3);
  auto b2 = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
  auto c2 = 2 * (q2 * q3 - q0 * q1);

  auto a3 = 2 * (q1 * q3 - q0 * q2);
  auto b3 = 2 * (q2 * q3 + q0 * q1);
  auto c3 = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

  // auto pos_dot = body_to_earth(u, w, a1, c1, a2, c2, a3, c3);
  const Vec3 pos_dot = Vec3(a1 * u + c1 * w, a2 * u + c2 * w, a3 * u + c3 * w);

  auto p = 0.0;
  auto qdot = 0;
  // vdot = -r * u + g * b3
  auto r = parms.g * b3 / u;

  auto udot = -q * w + a3 * parms.g + QS_m * Cx;
  auto wdot = q * u + c3 * parms.g + QS_m * Cz;
  // auto vel_dot = body_to_earth(udot, wdot);

  auto qmag = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
  auto lam = 1 - qmag;
  auto q0dot = -0.5 * (q1 * p + q2 * q + q3 * r) + lam * q0;
  auto q1dot = 0.5 * (q0 * p + q2 * r - q3 * q) + lam * q1;
  auto q2dot = 0.5 * (q0 * q - q1 * r + q3 * p) + lam * q2;
  auto q3dot = 0.5 * (q0 * r + q1 * q - q2 * p) + lam * q3;

  DerivState dstates;
  dstates[POS_X] = std::get<0>(pos_dot);
  dstates[POS_Y] = std::get<1>(pos_dot);
  dstates[POS_Z] = std::get<2>(pos_dot);
  dstates[VEL_U] = udot;
  dstates[VEL_W] = wdot;
  dstates[AVEL_Q] = qdot;
  dstates[QUATERNION] = q0dot;
  dstates[QUATERNION + 1] = q1dot;
  dstates[QUATERNION + 2] = q2dot;
  dstates[QUATERNION + 3] = q3dot;
  return dstates;
  // theta = sympy.asin(-a3 / qmag);
  // psi = sympy.acos(a1 / qmag / sympy.cos(theta)) * sympy.sign(a2);
  // phi = sympy.acos(c3 / qmag / sympy.cos(theta)) * sympy.sign(b3);
}

void system_model(State& state, double dt)
{
  DerivState deriv = system_ode(state);
  auto& [states, attitude] = state.data;
  for (int i = 0; i < 6; ++i)
  {
    states[i] += dt * deriv[i];
  }
  auto& quaternion = attitude.get_q().coeffs();
  double qs[4];
  for (int i = 0; i < 4; ++i)
  {
    qs[i] = quaternion[i] + dt * deriv[6 + i];
  }
  attitude =
      unscented::UnitQuaternion(Eigen::Quaterniond(qs[0], qs[1], qs[2], qs[3]));
}

//////////////////////////////////////////////////////////////////////////////
// Set up the measurement
//////////////////////////////////////////////////////////////////////////////

using Measurement = unscented::Compound<Pos, Pos, Pos, Vel>;

//////////////////////////////////////////////////////////////////////////////
// Set up the measurement model
//////////////////////////////////////////////////////////////////////////////

Measurement measurement_model(const State& state)
{
  auto& [states, attitude] = state.data;
  auto& u = states[VEL_U];
  auto& w = states[VEL_W];
  auto V_sq = u * u + w * w;
  auto V = sqrt(V_sq);
  return {states[POS_X], states[POS_Y], states[POS_Z], V};
}

int main()
{
}