#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <regex>
#include <string>
#include <tclap/CmdLine.h>
#include <valarray>
#include <vector>


using namespace std;


const string VERSION = "0.1.1";


struct Arguments {
  int n0;
  vector<double> birth_rate;
  double interaction_death_rate;
  vector<double> times;
};


class Cell {
public:
  Cell(vector<double> br, double q, float t_end)
    : birth_rates(br)
    , q(q)
    , t_end(t_end) {

    dt = t_end / (birth_rates.size() - 1);
  }

  // interpolate between a and b at coordinate x = [0, 1]
  double interpolate(double a, double b, double x) {
    return a*(1.0 - x) + b*x;
  }


  double get_birth_rate(double t) {
    if (t > t_end)
      return birth_rates.back();
    size_t i = t / dt;
    return interpolate(birth_rates[i], birth_rates[i+1], (t - i*dt)/dt);
  }

  vector<double> birth_rates;
  double dt;
  double q;
  double t_end;
};


template <typename TCell, typename TRng=std::mt19937>
class LB {
public:
  LB(TCell wt) {
    urd = std::uniform_real_distribution<double>(std::nextafter(0.0, 1.0), 1.0);
    std::random_device rd;
    rng.seed(rd());
    type_count = 1;
    X.resize(type_count);
    a.resize(type_count * 2);
    T.resize(type_count * 2);
    P.resize(type_count * 2);
    r.resize(type_count * 2);
    dt.resize(type_count * 2);
    X = 0;
    // a has to be calculated, so no reason to waste time initializing
    T = 0;
    P = 0;
    for (size_t i = 0; i < type_count; ++i) {
      cells.push_back(wt);
    }
  }

  void set_cell_count(size_t count) {
    X[0] = count;
  }


  auto get_cell_count() {
    return X[0];
  }

  auto get_time() {
    return t;
  }

  auto get_birth_rate(size_t cell_type) {
    return a[cell_type*2];
  }


  auto get_death_rate(size_t cell_type) {
    return a[cell_type*2 + 1];
  }



  void get_birth_rates() {
    for (size_t i = 0; i < type_count; ++i) {
      double rate = cells[i].get_birth_rate(t);
      double ai = X[i] * rate;
      a[i*2] = ai;
    }
  }
  void get_death_rates() {
    for (size_t i = 0; i < type_count; ++i) {
      double interaction = cells[i].q * cells[i].get_birth_rate(t); // Constant carrying capacity
      int sizemult = X.sum() - 1;
      a[i*2 + 1] = sizemult * X[i] * interaction;
    }
  }

  void update_rates() {
    get_death_rates(); // order is not important now (I think)!
    get_birth_rates();
  }

  void init() {
    for (auto &rr: r) {
      rr = urd(rng);
    }
    P = log(1.0/r);
    // calculate all propensity values (a)
    // These functions should probably be optimized with sfinae or something
    update_rates();
  }


  void simulate(double interval) {
    double t_end = t + interval;
    if (t == 0.0)
      init();

    while (t < t_end) {
      dt = (P - T) / a;

      size_t u = std::min_element(std::begin(dt), std::end(dt)) - std::begin(dt);
      double d = dt[u];
      t += d;
      int event_celltype = u/2;
      // std::cout << u << '\n';
      if (u%2 == 0) {
        ++X[event_celltype];
      } else {
        --X[event_celltype];
      }
      T = T + a*d;
      r[u] = urd(rng);
      P[u] += log(1.0/r[u]);
      update_rates();

    }
  }

private:
  double t = 0.0;
  size_t type_count;
  std::valarray<double> X;
  std::valarray<double> a;
  std::valarray<double> T;
  std::valarray<double> P;
  std::valarray<double> r;
  std::valarray<double> dt;

  std::vector<TCell> cells;

  TRng rng;
  std::uniform_real_distribution<double> urd;

};


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
    regex re_fp("\\d+\\.\\d+");
    while (regex_search(bstring, m_b, re_fp)) {
      for (auto x: m_b) {
        a.birth_rate.push_back(stod(x));
      }
      bstring = m_b.suffix().str();
    }
    a.interaction_death_rate = a_interaction_death_rate.getValue();
    string tstring = a_times.getValue();
    smatch m_t;
    while (regex_search(tstring, m_t, re_fp)) {
      for (auto x: m_t) {
        a.times.push_back(stod(x));
      }
      tstring = m_t.suffix().str();
    }

    // Sanity checks
    assert(a.n0 > 0);
    assert(a.interaction_death_rate >= 0.0);
    for (auto rate: a.birth_rate) {
      assert(rate >= 0.0);
    }
    for (size_t i = 1; i < a.times.size(); ++i) {
      assert(a.times[i-1] >= 0.0);
      // verify that they are sorted
      assert(a.times[i-1] <= a.times[i]);
    }


  } catch (TCLAP::ArgException &e) {
    cerr << "TCLAP Error: " << e.error() << endl << "\targ: " << e.argId() << endl;
    return 1;
  }

  // ### Simulation ### //

  Cell wt(a.birth_rate, a.interaction_death_rate, a.times.back());
  LB<Cell> lb(wt);
  lb.set_cell_count(a.n0);
  std::cout << "time\tsize\trate\n";
  double t_prev = 0.0;
  for (auto time: a.times) {
    lb.simulate(time - t_prev);
    t_prev = time;
    std::cout << time << '\t' << lb.get_cell_count() << '\t' << wt.get_birth_rate(time) << '\n';
  }
}

