#include "matplot/freestanding/axes_functions.h"
#include "unscented/primitives.h"
#include "unscented/ukf.h"

#include <matplot/matplot.h>
#include <iostream>
#include <random>
#include "DOP853.h"

using namespace tableau::integration;

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

using DerivState = std::vector<double>;

using Measurement = unscented::Compound<Pos, Pos, Pos, Vel>;

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

DerivState system_ode(const DerivState& state)
{
  auto& u = state[VEL_U];
  auto& w = state[VEL_W];
  auto& q = state[AVEL_Q];
  auto& q0 = state[QUATERNION + 0];
  auto& q1 = state[QUATERNION + 1];
  auto& q2 = state[QUATERNION + 2];
  auto& q3 = state[QUATERNION + 3];

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
  const auto qdot = 0;
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

  DerivState dstates = {std::get<0>(pos_dot),
                        std::get<1>(pos_dot),
                        std::get<2>(pos_dot),
                        udot,
                        wdot,
                        qdot,
                        q0dot,
                        q1dot,
                        q2dot,
                        q3dot};
  return dstates;
  // theta = sympy.asin(-a3 / qmag);
  // psi = sympy.acos(a1 / qmag / sympy.cos(theta)) * sympy.sign(a2);
  // phi = sympy.acos(c3 / qmag / sympy.cos(theta)) * sympy.sign(b3);
}

DerivState convert_state(const State& state)
{
  auto& [y, quaternion] = state.data;
  DerivState vstate;
  for (int i = 0; i < 6; ++i)
  {
    vstate.push_back(y[i]);
  }
  auto& q = quaternion.get_q();
  vstate.push_back(q.w());
  vstate.push_back(q.x());
  vstate.push_back(q.y());
  vstate.push_back(q.z());
  return vstate;
}

void set_state(State& state, const DerivState& vstate)
{
  auto& [y, quaternion] = state.data;
  for (int i = 0; i < 6; ++i)
  {
    y[i] = vstate[i];
  }
  quaternion = unscented::UnitQuaternion(
      Eigen::Quaterniond(vstate[QUATERNION + 0], vstate[QUATERNION + 1],
                         vstate[QUATERNION + 2], vstate[QUATERNION + 3]));
}

void write(const DerivState& x)
{
  for (int i = 0; i < 10; ++i)
  {
    std::cout << x[i] << ", ";
  }
  std::cout << std::endl;
}

void system_model(State& state, double dt)
{
  double t0 = 0.0;
  DerivState y0 = convert_state(state);
  DerivState dy0 = system_ode(y0);
  write(y0);

  DOP853Config<DerivState> cfg;
  cfg.derivative = [](const DerivState& y, double) { return system_ode(y); };
  DOP853Integrator<DerivState> integrator(cfg);
  auto result = integrator.integrate(t0, y0, dt,
                                     DOP853Tolerance::scalar(1.0e-12, 1.0e-12));
  set_state(state, result.y);
}

//////////////////////////////////////////////////////////////////////////////
// Set up the measurement
//////////////////////////////////////////////////////////////////////////////

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
  // Initialize the UKF and set the weights
  using UKF = unscented::UKF<State, Measurement>;
  UKF ukf;
  ukf.set_weight_coefficients(0.1, 2.0, -1.0);

  // Simulation parameters
  const auto DT = 1.0; // seconds

  UKF::N_by_N Q;
  Q.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01;
  ukf.set_process_covariance(Q);

  UKF::M_by_M R;
  R.diagonal() << 4.0, 4.0, 4.0, 1.0;
  ukf.set_measurement_covariance(R);

  // Set initial state estimate and its covariance
  State true_state({0, 0, 0, 20, 0, 0}, unscented::UnitQuaternion());
  //   State initial_state_estimate(0, 90, 1100, 0);
  State initial_state_estimate({0, 0, 0, 20, 0, 0},
                               unscented::UnitQuaternion());
  ukf.set_state(initial_state_estimate);
  UKF::N_by_N P;
  P.diagonal() << 10.0, 10.0, 10.0, 1.0, 1.0, 0.1, 1.0, 1.0, 1.0;
  ukf.set_state_covariance(P);

  double t = 0.0;

  while (t < 200.0)
  {
    system_model(true_state, DT);
    t += DT;
  }
}