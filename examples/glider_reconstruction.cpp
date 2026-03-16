#include "matplot/freestanding/axes_functions.h"
#include "unscented/primitives.h"
#include "unscented/ukf.h"

#include <matplot/matplot.h>
#include <iostream>
#include <random>
#include "DOP853.h"

using namespace tableau::integration;

#include <boost/json.hpp>
#include <fstream>
#include <string>
namespace json = boost::json;

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

const GliderAero parms;

DerivState system_ode(const DerivState& state)
{
  auto& u = state[VEL_U];
  auto& w = state[VEL_W];
  auto& q = state[AVEL_Q];
  auto& q0 = state[QUATERNION + 0];
  auto& q1 = state[QUATERNION + 1];
  auto& q2 = state[QUATERNION + 2];
  auto& q3 = state[QUATERNION + 3];

  const auto V_sq = u * u + w * w;
  const auto V = sqrt(V_sq);
  const auto alpha = atan(w / u);

  auto CL = parms.dcldalpha * alpha + parms.CL0;
  // CL = sympy.Min(1.5, sympy.Max(0.0, cl_raw))
  auto CD = parms.CD0 + parms.k * CL * CL;
  // CLL = sympy.Min(CL, 1.0)
  auto CLL = CL;

  const auto Q = 0.5 * parms.rho * V_sq;
  const auto QS_m = Q * parms.S / parms.m;
  const auto sa = w / V; // sin(alpha)
  const auto ca = u / V; // cos(alpha)
  const auto Cx = CLL * sa - CD * ca;
  const auto Cz = -CLL * ca - CD * sa;

  // orientation
  const Eigen::Matrix3d R =
      Eigen::Quaterniond(q0, q1, q2, q3).toRotationMatrix();
  const Eigen::Vector3d U(u, 0, w);
  const Eigen::Vector3d pos_dot = R * U;

  const auto p = 0.0;
  const auto qdot = 0;
  // vdot = -r * u + g * b3 = 0
  const auto r = parms.g * R(2, 1) / u;

  auto udot = -q * w + R(2, 0) * parms.g + QS_m * Cx;
  auto wdot = q * u + R(2, 2) * parms.g + QS_m * Cz;

  auto qmag = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
  auto lam = 1 - qmag;
  auto q0dot = -0.5 * (q1 * p + q2 * q + q3 * r) + lam * q0;
  auto q1dot = 0.5 * (q0 * p + q2 * r - q3 * q) + lam * q1;
  auto q2dot = 0.5 * (q0 * q - q1 * r + q3 * p) + lam * q2;
  auto q3dot = 0.5 * (q0 * r + q1 * q - q2 * p) + lam * q3;

  return DerivState({pos_dot(0), pos_dot(1), pos_dot(2), udot, wdot, qdot,
                     q0dot, q1dot, q2dot, q3dot});
}

const DerivState convert_state(const State& state)
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

template <typename T>
int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

const Eigen::Vector3d get_euler(const State& state)
{
  auto& [y, quaternion] = state.data;
  const Eigen::Matrix3d R = quaternion.get_q().toRotationMatrix();
  double theta = asin(-R(2, 0));
  double psi = acos(R(0, 0) / cos(theta)) * sign(R(1, 0));
  double phi = acos(R(2, 2) / cos(theta)) * sign(R(2, 1));
  return Eigen::Vector3d(phi, theta, psi);
  //  * 180.0f / M_PI;
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

void write(const State& state)
{
  auto& x = convert_state(state);
  for (int i = 0; i < 6; ++i)
  {
    std::cout << x[i] << ", ";
  }
  auto& euler = get_euler(state);
  for (int i = 0; i < 3; ++i)
  {
    std::cout << euler[i] << ", ";
  }
  std::cout << std::endl;
}

void system_model(State& state, double dt)
{
  double t0 = 0.0;
  DerivState y0 = convert_state(state);
  DerivState dy0 = system_ode(y0);

  DOP853Config<DerivState> cfg;
  cfg.derivative = [](const DerivState& y, double) { return system_ode(y); };
  DOP853Integrator<DerivState> integrator(cfg);
  auto result = integrator.integrate(t0, y0, dt,
                                     DOP853Tolerance::scalar(1.0e-12, 1.0e-12));
  if (result.y[VEL_U] < 1.0)
  {
    result.y[VEL_U] = 1.0;
  }
  set_state(state, result.y);
}

/// @brief /////////////////////////////////////////////////////////

// Utility function to read file content into a string
std::string read_file(const char* filename)
{
  std::ifstream is(filename);
  if (!is)
  {
    throw std::runtime_error(std::string("Cannot open file: ") + filename);
  }
  std::string s;
  is.seekg(0, std::ios::end);
  s.resize(is.tellg());
  is.seekg(0, std::ios::beg);
  is.read(&s[0], s.size());
  return s;
}

State get_initial_state_estimate(const double x, const double y, const double z,
                                 const double U, const double hdg,
                                 const double pitch, const double bank)
{
  Eigen::Quaterniond init_quaternion =
      Eigen::AngleAxisd(hdg, Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(bank, Eigen::Vector3d::UnitX());
  return State({x, y, z, U, 0.0, 0.0},
               unscented::UnitQuaternion(init_quaternion));
}

std::vector<Measurement> load_encounter(State& initial_state_estimate)
{
  std::vector<Measurement> measurements;
  try
  {
    // Read the file into a string
    auto const json_data = read_file(
        "/home/jmw/Desktop/Projects/FlightReconstruction/encounter_00000.json");

    // Parse the JSON string into a value
    json::value jv = json::parse(json_data);
    auto trace = jv.get_object().at("aircraft").at(0).at("trace").as_array();
    int count = 0;
    for (auto& line : trace)
    {
      if (count == 2)
      {
        auto& [x, y, z, U] = measurements[0].data;
        const double hdg = line.at("hdg").as_double() * M_PI / 180.0;
        const double pitch = line.at("pitch").as_double() * M_PI / 180.0;
        const double bank = line.at("bank").as_double() * M_PI / 180.0;
        initial_state_estimate = get_initial_state_estimate(
            x.value, y.value, z.value, U.value, hdg, pitch, bank);
      }
      if (count > 0)
      {
        Measurement measurement;
        auto& [y, x, z, U] = measurement.data;
        x.value = line.at("x").as_double();
        y.value = line.at("y").as_double();
        z.value = -line.at("alt_gps").as_double();
        U.value = line.at("v_tas").as_double();
        measurements.push_back(measurement);
      }
      count++;
    }
  }
  catch (std::exception const& e)
  {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
  // json::value jv = json::parse(json_data);
  return measurements;
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

using UKF = unscented::UKF<State, Measurement>;

UKF get_ukf(const State& initial_state_estimate, const double DT)
{
  UKF ukf;
  ukf.set_weight_coefficients(0.1, 2.0, -1.0);

  UKF::N_by_N Q;
  Q.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01;
  ukf.set_process_covariance(Q);

  UKF::M_by_M R;
  R.diagonal() << 4.0, 4.0, 4.0, 1.0;
  ukf.set_measurement_covariance(R);

  ukf.set_state(initial_state_estimate);
  UKF::N_by_N P;
  P.diagonal() << 10.0, 10.0, 10.0, 1.0, 1.0, 0.1, 1.0, 1.0, 1.0;
  ukf.set_state_covariance(P);

  return ukf;
}

int main()
{
  // Simulation parameters
  const auto DT = 1.0; // seconds

  State initial_state_estimate;
  std::vector<Measurement> measurements =
      load_encounter(initial_state_estimate);

  auto ukf = get_ukf(initial_state_estimate, DT);
  double t = 0.0;

  for (auto& meas : measurements)
  {
    // system_model(true_state, DT);
    // auto meas = measurement_model(true_state);
    // auto& [range, elevation] = meas.data;
    // range.value += range_noise(gen);
    // elevation = unscented::Angle(elevation.get_angle() +
    // elevation_noise(gen));

    // Update the filter estimates
    ukf.predict(system_model, DT);
    ukf.correct(measurement_model, meas);
    ukf.smooth();
    t += DT;
    const auto& est_state = ukf.get_state();
    write(est_state);
  }
}