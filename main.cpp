#include <rubikscube.hpp>
#include <random>
#include <chrono>

random_device rd;
mt19937 mt(rd());

int get_random(int low, int high) {
  uniform_int_distribution<int> range(low, high);
  return range(mt);
}

void scramble(RubiksCube& cube) {
  const int NUM_MOVES = 40;
  const int POSSIBLE_MOVES = 18;
  for (int i = 0; i < NUM_MOVES; ++i) {
    int move = get_random(0, POSSIBLE_MOVES - 1);
    cube.apply_move(static_cast<Move>(move));
  }
}

double time_for_solving() {
  RubiksCube cube;
  scramble(cube);
  const auto t_beg = std::chrono::steady_clock::now();
  KociembaSolver solver;
  solver.solve(cube);
  const auto t_end = std::chrono::steady_clock::now();
  const std::chrono::duration<double> diff = t_end - t_beg;
  return diff.count();
}

int main() {
  const int TESTS = 1000;
  double max_time = 0;
  double min_time = 100;
  double average_time = 0;
  for (int i = 0; i < TESTS; ++i) {
    double time = time_for_solving();
    max_time = max(max_time, time);
    min_time = min(min_time, time);
    average_time += time;
  }
  average_time /= TESTS;
  cout << "Count of tests = " << TESTS << '\n';
  cout << "Average solving time = " << average_time << '\n';
  cout << "Minimum solving time = " << min_time << '\n';
  cout << "Maximum solving time = " << max_time << '\n';
}