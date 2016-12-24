#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <queue>
#include <chrono>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <cstring>
#include <utility>
#include <ctime>

namespace
{
  typedef std::vector<std::pair<int, double>> Distribution;
  typedef std::vector<Distribution> RectDistribution;

  void generate_distribution(size_t n, bool desc, Distribution& distr) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine gen(seed);
    std::uniform_int_distribution<int> rnd(1, n << 2);

    int sum = 0;
    std::vector<int> vals(n);
    for (int i = 0; i < n; ++i) {
      vals[i] = rnd(gen);
      sum += vals[i];
    }

    if (desc) {
      std::sort(vals.rbegin(), vals.rend());
    } else {
      std::sort(vals.begin(), vals.end());
    }

    for (int i = 0; i < n; ++i) {
      distr[i] = std::make_pair(i, 1.0 * vals[i] / sum);
    }
  }

  void sample_from_distribution(const Distribution& distr, size_t n, size_t sample_size, std::ostream& os, std::ostream& log) {
    log << "sample_from_distribution: " << "\n";

    double* prefix_sum = new double[n];

    prefix_sum[0] = distr[0].second;
    for (int i = 1; i < n - 1; ++i) {
      prefix_sum[i] = prefix_sum[i - 1] + distr[i].second;
    }
    prefix_sum[n - 1] = 1.0;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine gen(seed);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);

    os << sample_size << "\n";

    clock_t start = std::clock();
    for (int i = 0; i < sample_size; ++i) {
      double val = rnd(gen);
      for (int j = 0; j < n; ++j) {
        if (val <= prefix_sum[j]) {
          os << distr[j].first << "\n";
          break;
        }
      }
    }
    log << "-time elapsed: " << double(std::clock() - start) / CLOCKS_PER_SEC << "s\n";

    delete[] prefix_sum;
  }

  void sample_from_distribution_bin(const Distribution& distr, size_t n, size_t sample_size, std::ostream& os, std::ostream& log) {
    log << "sample_from_distribution_bin: " << "\n";

    double* prefix_sum = new double[n];

    prefix_sum[0] = distr[0].second;
    for (int i = 1; i < n - 1; ++i) {
      prefix_sum[i] = prefix_sum[i - 1] + distr[i].second;
    }
    prefix_sum[n - 1] = 1.0;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine gen(seed);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);

    os << sample_size << "\n";

    clock_t start = std::clock();
    for (int i = 0; i < sample_size; ++i) {
      double val = rnd(gen);
     
      int pos = -1; 

      int low = 0;
      int high = n - 1;
      while (low <= high) {
        int mid = (low + high) / 2;
        if (val < prefix_sum[mid]) {
          high = mid - 1;
        } else if (val > prefix_sum[mid]){
          low = mid + 1;
        } else {
          pos = mid;
          break;
        }
      }
      if (pos == -1) {
        pos = low;
      }

      os << distr[pos].first << "\n";
    }
    log << "-time elapsed: " << double(std::clock() - start) / CLOCKS_PER_SEC << "s\n";

    delete[] prefix_sum;
  }

  bool func_cmp(const std::pair<int, double>& left, const std::pair<int, double>& right) {
    return left.second < right.second;
  }

  void sample_from_distribution_rect(const Distribution& distr, size_t n, size_t sample_size, std::ostream& os, std::ostream& log) {
    log << "sample_from_distribution_rect: " << "\n";

    Distribution _distr;
    _distr.assign(distr.begin(), distr.end());
    std::sort(_distr.begin(), _distr.end(), func_cmp);

    RectDistribution rect_distr(n, Distribution(2));

    for (int i = 0; i < n; ++i) 
    {
      rect_distr[i][0] = _distr[i];
      rect_distr[i][1] = std::make_pair(_distr[n - 1].first, 1.0 / n - _distr[i].second);

      if (i == n - 1) {
        break;
      }

      double diff = _distr[n - 1].second - rect_distr[i][1].second;

      int pos = -1;

      int low = i + 1;
      int high = n - 2;
      while (low <= high) {
        int mid = (low + high) / 2;
        if (diff < _distr[mid].second) {
          high = mid - 1;
        } else if (diff > _distr[mid].second){
          low = mid + 1;
        } else {
          pos = mid;
          break;
        }
      }
      if (pos == -1) {
        pos = low;
      }

      for (int j = n - 1; j > pos; --j) {
        _distr[j] = _distr[j - 1];
      }
      _distr[pos] = std::make_pair(rect_distr[i][1].first, diff);
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine gen(seed);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    std::uniform_real_distribution<double> rnd2(0.0, 1.0 / n);

    os << sample_size << "\n";

    clock_t start = std::clock();
    for (int i = 0; i < sample_size; ++i) {
      double val = rnd(gen);
      double val2 = rnd2(gen);

      int id_x = (n - 1) * val;
      int id_y = (val2 <= rect_distr[id_x][0].second ? 0 : 1);

      os << rect_distr[id_x][id_y].first << "\n";
    }
    log << "-time elapsed: " << double(std::clock() - start) / CLOCKS_PER_SEC << "s\n";
  }

  void run_test(const Distribution& distr, size_t n, size_t sample_size, const char* of) {
    std::ofstream log;
    log.open("distr.log", std::ofstream::out | std::ofstream::app);
    log << "\n\n\nTest configuration: \n" << "-n: " << n << "\n-sample_size: " << sample_size << "\n-output: " << of << "\n\n";

    std::ofstream ofs;
    ofs.open(of, std::ofstream::out);
    sample_from_distribution(distr, n, sample_size, ofs, log);
    sample_from_distribution_bin(distr, n, sample_size, ofs, log);
    sample_from_distribution_rect(distr, n, sample_size, ofs, log);

    ofs.close();
    log.close();
  }
}

int main(int argc, char** argv)
{
  Distribution distr100(100);
  generate_distribution(100, false, distr100);
  run_test(distr100, 100, 1000000, "distr100.dat");

  Distribution distr400(400);
  generate_distribution(400, true, distr400);
  run_test(distr400, 400, 1000000, "distr400.dat");

  Distribution distr700(700);
  generate_distribution(700, false, distr700);
  run_test(distr700, 700, 1000000, "distr700.dat");

  Distribution distr1000(1000);
  generate_distribution(1000, true, distr1000);
  run_test(distr1000, 1000, 1000000, "distr1000.dat");

  return 0;
}