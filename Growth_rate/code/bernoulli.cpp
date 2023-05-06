#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <regex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <tclap/CmdLine.h>


using namespace std;


const string VERSION = "0.1.0";


struct Arguments {
  int n0;
  vector<double> birth_rate;
  double interaction_death_rate;
  vector<double> times;
};


static vector<double> birth_rate;
static double t_end;
static double interaction;


// interpolate between a and b at coordinate x = [0, 1]
double interpolate(double a, double b, double x) {
  return a*(1.0 - x) + b*x;
}


double birth_rate_function(double t) {
  double dt = t_end / (birth_rate.size() - 1);
  size_t i = t / dt;
  return interpolate(birth_rate[i], birth_rate[i+1], (t - i*dt)/dt);
}


double birth_rate_integral(double a, double b) {

  bool negate = false;
  if (a > b) {
    negate = true;
    swap(a, b);
  }

  double integral = 0.0;
  double dt = t_end / (birth_rate.size() - 1);
  size_t start_segment = a / dt;
  size_t end_segment = b / dt;

  // If b is equal to end_time, it is still in the final segment
  if (end_segment == birth_rate.size() - 1)
    --end_segment;

  double low_birthrate = interpolate(birth_rate[start_segment], birth_rate[start_segment + 1], (a - dt*start_segment)/dt);
  double high_birthrate = interpolate(birth_rate[end_segment], birth_rate[end_segment + 1], (b - dt*end_segment)/dt);

  if (start_segment == end_segment) {
    integral += low_birthrate * (b - a);
    integral += (high_birthrate - low_birthrate) * (b - a) / 2.0;
  } else {
    // first and final segment
    integral += low_birthrate * (dt*(start_segment + 1) - a);
    integral += (birth_rate[start_segment + 1] - low_birthrate) * (dt*(start_segment + 1) - a) / 2.0;
    integral += birth_rate[end_segment] * (b - dt*end_segment);
    integral += (high_birthrate - birth_rate[end_segment]) * (b - dt*end_segment) / 2.0;

    // intervening segments
    for (size_t i = start_segment + 1; i < end_segment; ++i) {
      integral += birth_rate[i] * dt;
      integral += (birth_rate[i+1] - birth_rate[i]) * dt / 2.0;
    }
  }

  if (negate)
    return -integral;

  return integral;
}



double denominator_function(double t, [[maybe_unused]] void *params) {
  return -interaction * birth_rate_function(t) * exp(birth_rate_integral(0.0, t));
}


int main(int argc, char **argv) {

  // ### Argument parsing ### //

  Arguments a;
  try {
    TCLAP::CmdLine cmd("General treatment simulator", ' ', VERSION);

    TCLAP::ValueArg<int> a_n0("n", "n0", "Starting cell count", true, 100, "integer", cmd);
    TCLAP::ValueArg<string> a_birth_rate("b", "birth-rate", "Birth rate", true, "", "[0, 1, 2, ...]", cmd);
    TCLAP::ValueArg<double> a_interaction_death_rate("q", "interaction-death_rate", "Interaction Death rate", true, 100, "double", cmd);
    TCLAP::ValueArg<string> a_times("t", "measure-times", "Times to measure population size", true, "", "[0, 1, 2, ...]", cmd);

    cmd.parse(argc, argv);

    a.n0 = a_n0.getValue();
    string bstring = a_birth_rate.getValue();
    smatch m_b;
    regex re_b("[-]?\\d+\\.\\d+");
    while (regex_search(bstring, m_b, re_b)) {
      for (auto x: m_b) {
        a.birth_rate.push_back(stod(x));
      }
      bstring = m_b.suffix().str();
    }
    if (a.birth_rate.size() == 1) {
      // we always need at least 2 values to define the piecewise linear curve
      // interpret a single value as a constant line
      a.birth_rate.push_back(a.birth_rate[0]);
    }
    a.interaction_death_rate = a_interaction_death_rate.getValue();
    string tstring = a_times.getValue();
    smatch m_t;
    while (regex_search(tstring, m_t, re_b)) {
      for (auto x: m_t) {
        a.times.push_back(stod(x));
      }
      tstring = m_t.suffix().str();
    }

    // Sanity checks
    assert(a.n0 > 0);
    for (size_t i = 1; i < a.times.size(); ++i) {
      assert(a.times[i-1] >= 0.0);
      // verify that they are sorted
      assert(a.times[i-1] <= a.times[i]);
    }
    assert(a.interaction_death_rate >= 0);

  } catch (TCLAP::ArgException &e) {
    cerr << "TCLAP Error: " << e.error() << endl << "\targ: " << e.argId() << endl;
    return 1;
  }

  // ### Calculation ### //

  // setup globals
  t_end = a.times.back();
  interaction = a.interaction_death_rate;
  birth_rate = a.birth_rate;

  // Logistic growth with a variable birthrate is a bernoulli differential equation with the following solution
  // f(t) = e^( integral_0^t a(ξ) dξ)/(c_1 - integral_0^t-b e^( integral_0^ζ a(ξ) dξ) dζ)
  // where a is the birthrate and b is the interaction factor

  gsl_integration_workspace *workspace;

  workspace = gsl_integration_workspace_alloc(1000);

  gsl_function denominator_integral;
  denominator_integral.function = &denominator_function;
  denominator_integral.params = nullptr;

  std::cout << "time\tsize\trate\n";
  cout.precision(numeric_limits<double>::max_digits10);

  // it is possible for the integration to fail on numerical errors
  // in that case, we want the program to keep running while raising the error limits
  gsl_set_error_handler_off();

  for (auto time: a.times) {
    double t = time;
    // First, find the numerator integral_0^t a(ξ) dξ
    double numerator = birth_rate_integral(0.0, t);

    // calculate integral in denominator with gsl
    double d_int;
    double d_int_err;

    int status;
    double tolerance = 1e-7;
    do {
      status = gsl_integration_qag(&denominator_integral, 0.0, t, 0, tolerance, 1000, 6, workspace, &d_int, &d_int_err);
      if (status) {
        tolerance *= 10.0;
      }
    } while (status != 0);

    // Finally, find population size at time t
    double N = exp(numerator) / (1.0/a.n0 - d_int);

    cout << t << '\t' << N << '\t' << birth_rate_function(t) << endl;
  }

  // reenable in case the workspace free fails
  gsl_set_error_handler(NULL);

  gsl_integration_workspace_free(workspace);
}
