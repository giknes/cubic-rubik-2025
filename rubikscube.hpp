#include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <assert.h>
#include <cstring>

using namespace std;

// ================ combinatorial functions ================
class Combinatorics {
private:
  static const int N = 12;
  static const int K = 4;
  array<array<int, K + 1>, N + 1> comb;

public:
  static Combinatorics& getInstance() {
    static Combinatorics instance;
    return instance;
  };
  Combinatorics() {
    // calculate binomial coefficient  C(n,k)
    for (int n = 0; n <= N; n++) {
      for (int k = 0; k <= K; k++) {
        if (k == 0 || k == n) comb[n][k] = 1;
        else if (k > n) comb[n][k] = 0;
        else comb[n][k] = comb[n - 1][k - 1] + comb[n - 1][k];
      }
    }
  }

  // encode position of 4 middle slices
  int getCombinationIndex(const array<uint8_t, 12>& edges) const {
    int index = 0;
    int k = K;
    for (int i = N - 1; i >= 0 && k > 0; i--) {
      if (edges[i] >= 8) {
        index += comb[i][k];
        k--;
      }
    }
    return index;
  }

  // encode permutation using factorial number system
  // index = (p[i] - 1)! * (count of elements to the right and smaller that p[i])
  template <size_t N>
  uint32_t get_permutation_index(const std::array<uint8_t, N>& elements) {
    uint32_t index = 0;
    for (size_t i = 0; i < N; ++i) {
      int count = 0;
      for (size_t j = i + 1; j < N; ++j) {
        if (elements[j] < elements[i]) count++;
      }
      index = index * (N - i) + count;
    }
    return index;
  }
};


// All possible moves of faces
enum Move {
  U1, U2, U3, D1, D2, D3,
  R1, R2, R3, L1, L2, L3,
  F1, F2, F3, B1, B2, B3,
  NO_MOVE
};

const vector<string> move_names = {
    "U", "U2", "U'", "D", "D2", "D'",
    "R", "R2", "R'", "L", "L2", "L'",
    "F", "F2", "F'", "B", "B2", "B'"
};

// =============================================
// Rubik's Cube class
// =============================================

class RubiksCube {
private:
  // indexation of edges and corners
  enum EdgePositions {
    UR, UF, UL, UB, DR, DF, DL, DB, FR, FL, BL, BR
  };

  enum CornerPositions {
    URF, UFL, ULB, UBR, DFR, DLF, DBL, DRB
  };

public:
  array<uint8_t, 12> edges;
  array<uint8_t, 12> edge_orientations; // 0 - not flipped, 1 - flipped
  array<uint8_t, 8> corners;
  array<uint8_t, 8> corner_orientations; // 0 - not twisted, 1 - clockwise, 2 - counter-clockwise
  RubiksCube() {
    // Solver state of cube
    for (uint8_t i = 0; i < 12; i++) {
      edges[i] = i;
      edge_orientations[i] = 0;
    }
    for (uint8_t i = 0; i < 8; i++) {
      corners[i] = i;
      corner_orientations[i] = 0;
    }
  }

  // applying a move to a cube
  void apply_move(Move move) {
    switch (move) {
      // Up face (does not affect on orientation)
    case U1:
      rotate_edges(0, 1, 2, 3); rotate_corners(0, 1, 2, 3); break;
    case U2:
      swap(edges[0], edges[2]);
      swap(edges[1], edges[3]);
      swap(edge_orientations[0], edge_orientations[2]);
      swap(edge_orientations[1], edge_orientations[3]); 
      swap(corners[0], corners[2]);
      swap(corners[1], corners[3]);
      swap(corner_orientations[0], corner_orientations[2]);
      swap(corner_orientations[1], corner_orientations[3]); 
      break;
    case U3: rotate_edges(3, 2, 1, 0); rotate_corners(3, 2, 1, 0); break;

      // Down face (does not affect on orientation)
    case D1:   rotate_edges(7, 6, 5, 4); rotate_corners(7, 6, 5, 4); ; break;
    case D2: 
      swap(edges[4], edges[6]);
      swap(edges[5], edges[7]);
      swap(edge_orientations[4], edge_orientations[6]);  
      swap(edge_orientations[5], edge_orientations[7]);  
      swap(corners[4], corners[6]);
      swap(corners[5], corners[7]);
      swap(corner_orientations[4], corner_orientations[6]);  
      swap(corner_orientations[5], corner_orientations[7]);  
      break;
    case D3:     rotate_edges(4, 5, 6, 7); rotate_corners(4, 5, 6, 7); ; break;

      // Right face (affects on orientation: edges flip, corners which change UD slice +2, otherwise +1)
    case R1:
      rotate_edges(0, 11, 4, 8, { 1, 1, 1, 1 });
      rotate_corners(0, 3, 7, 4, { 2, 1, 2, 1 });
      break;
    case R2:
      swap(edges[0], edges[4]);
      swap(edges[8], edges[11]);
      swap(edge_orientations[0], edge_orientations[4]);  
      swap(edge_orientations[8], edge_orientations[11]);  
      swap(corners[0], corners[7]);
      swap(corners[3], corners[4]);
      swap(corner_orientations[0], corner_orientations[7]); 
      swap(corner_orientations[3], corner_orientations[4]);  
      break;
    case R3:
      rotate_edges(8, 4, 11, 0, { 1, 1, 1, 1 });
      rotate_corners(4, 7, 3, 0, { 1, 2, 1, 2 });
      break;

      // Left face (affects on orientation: edges flip, corners which change UD slice +2, otherwise +1)
    case L1:
      rotate_edges(2, 9, 6, 10, { 1, 1, 1, 1 });
      rotate_corners(2, 1, 5, 6, { 2, 1, 2, 1 });
      break;
    case L2:
      swap(edges[2], edges[6]);
      swap(edges[9], edges[10]);
      swap(edge_orientations[2], edge_orientations[6]);
      swap(edge_orientations[9], edge_orientations[10]);
      swap(corners[2], corners[5]);
      swap(corners[1], corners[6]);
      swap(corner_orientations[2], corner_orientations[5]);
      swap(corner_orientations[1], corner_orientations[6]);
      break;
      break;
    case L3:
      rotate_edges(10, 6, 9, 2, { 1, 1, 1, 1 });
      rotate_corners(6, 5, 1, 2, { 1, 2, 1, 2 });
      break;

      // Front face (affects on orientation: corners which change UD slice +2, otherwise +1)
    case F1:
      rotate_edges(1, 8, 5, 9, { 0, 0, 0, 0 });
      rotate_corners(1, 0, 4, 5, { 2, 1, 2, 1 });
      break;
    case F2:
      swap(edges[1], edges[5]);
      swap(edges[8], edges[9]);
      swap(edge_orientations[1], edge_orientations[5]);  
      swap(edge_orientations[8], edge_orientations[9]);
      swap(corners[0], corners[5]);
      swap(corners[1], corners[4]);
      swap(corner_orientations[0], corner_orientations[5]); 
      swap(corner_orientations[1], corner_orientations[4]); 
      break;
    case F3:
      rotate_edges(9, 5, 8, 1, { 0, 0, 0, 0 });
      rotate_corners(5, 4, 0, 1, { 1, 2, 1, 2 });
      break;

      // Back face (affects on orientation: corners which change UD slice +2, otherwise +1)
    case B1:
      rotate_edges(3, 10, 7, 11, { 0, 0, 0, 0 });
      rotate_corners(2, 6, 7, 3, { 1, 2, 1, 2 });
      break;
    case B2:
      swap(edges[3], edges[7]);   
      swap(edges[10], edges[11]); 
      swap(edge_orientations[3], edge_orientations[7]);    
      swap(edge_orientations[10], edge_orientations[11]); 
      swap(corners[2], corners[7]);  
      swap(corners[3], corners[6]); 
      swap(corner_orientations[2], corner_orientations[7]); 
      swap(corner_orientations[3], corner_orientations[6]);
      break;
    case B3:
      rotate_edges(11, 7, 10, 3, { 0, 0, 0, 0 });
      rotate_corners(3, 7, 6, 2, { 2, 1, 2, 1 });
      break;

    default: break;
    }
  }

private:
  // Helper functions for rotations
  template<class T>
  void rotate(T& a, T& b, T& c, T& d) {
    T temp = d;
    d = c;
    c = b;
    b = a;
    a = temp;
  }

  void rotate_edges(int a, int b, int c, int d, const vector<uint8_t>& flips = {}) {
    rotate(edges[a], edges[b], edges[c], edges[d]);
    rotate(edge_orientations[a], edge_orientations[b], edge_orientations[c], edge_orientations[d]);

    if (!flips.empty()) {
      edge_orientations[a] = (edge_orientations[a] + flips[0]) % 2;
      edge_orientations[b] = (edge_orientations[b] + flips[1]) % 2;
      edge_orientations[c] = (edge_orientations[c] + flips[2]) % 2;
      edge_orientations[d] = (edge_orientations[d] + flips[3]) % 2;
    }
  }

  void rotate_corners(int a, int b, int c, int d, const vector<uint8_t>& twists = {}) {
    rotate(corners[a], corners[b], corners[c], corners[d]);
    rotate(corner_orientations[a], corner_orientations[b], corner_orientations[c], corner_orientations[d]);

    if (!twists.empty()) {
      corner_orientations[a] = (corner_orientations[a] + twists[0]) % 3;
      corner_orientations[b] = (corner_orientations[b] + twists[1]) % 3;
      corner_orientations[c] = (corner_orientations[c] + twists[2]) % 3;
      corner_orientations[d] = (corner_orientations[d] + twists[3]) % 3;
    }
  }

public:
  bool is_solved() const {
    for (uint8_t i = 0; i < 12; i++)
      if (edges[i] != i || edge_orientations[i] != 0) return false;
    for (uint8_t i = 0; i < 8; i++)
      if (corners[i] != i || corner_orientations[i] != 0) return false;
    return true;
  }

  // Check for phase 1 (G1 <= H)
  bool in_phase1() const {
    // Check for orientation
    for (uint8_t i = 0; i < 12; i++)
      if (edge_orientations[i] != 0) return false;
    for (uint8_t i = 0; i < 8; i++)
      if (corner_orientations[i] != 0) return false;

    // Check for middle-edges (they must locate in middle slice)
    for (uint8_t i = 8; i < 12; i++)
      if (edges[i] < 8 || edges[i] > 11) return false;

    return true;
  }

  // ================ Phase 1 hash functions ================

  // 1. Edge orientation hash (12 bit, 0-4095)
  uint32_t get_edge_orientation_hash() const {
    uint32_t hash = 0;
    for (int i = 0; i < 12; ++i) {
      hash |= (1 << i) * edge_orientations[i];
    }
    return hash;
  }

  // 2. Corner orientation hash (7 elements define last one due to parity, 0-6560)
  uint32_t get_corner_orientation_hash() const {
    uint32_t hash = 0;
    for (int i = 0; i < 7; ++i) {
      hash = hash * 3 + corner_orientations[i];
    }
    return hash;
  }

  // 3. Middle-slice edges (0-494)
  uint32_t get_middle_comb_hash() const {
    return Combinatorics::getInstance().getCombinationIndex(edges);
  }

  // ================ Phase 2 hash functions ================

  // 4. Corner position hash (0-40319)
  uint32_t get_corner_position_hash() const {
    return Combinatorics::getInstance().get_permutation_index(corners);
  }

  // 5. UD-edges position hash (0-40319)
  uint32_t get_ud_edges_position_hash() const {
    array<uint8_t, 8> reduced_ud_edges;
    copy(edges.begin(), edges.begin() + 8, reduced_ud_edges.begin());
    return Combinatorics::getInstance().get_permutation_index(reduced_ud_edges);
  }

  // 6. Middle edges postion hash (0-23)
  uint32_t get_middle_edges_position_hash() const {
    array<uint8_t, 4> reduced_middle_edges;
    copy(edges.begin() + 8, edges.begin() + 12, reduced_middle_edges.begin());
    return Combinatorics::getInstance().get_permutation_index(reduced_middle_edges);
  }
};


// ====================================================
// Class for implementation of the Kociemba's algorithm
// ====================================================

class KociembaSolver {
private:
  static constexpr uint8_t INIT_VAL = 200;
  static const uint32_t EDGE_ORIENT = 4096;
  static const uint32_t CORNER_ORIENT = 6561;
  static const uint32_t MIDDLE_EDGES = 495;
  static const uint32_t UD_EDGES_POS = 40320;
  static const uint32_t CORNERS_POS = 40320;
  static const uint32_t MIDDLE_EDGES_POS = 24;
  // Pruning tables
  // Phase 1
  array<uint8_t, EDGE_ORIENT> edges_ori_pruning;
  array<uint8_t, CORNER_ORIENT> corners_ori_pruning;
  array<uint8_t, MIDDLE_EDGES> middle_edges_comb_pruning;

  // Phase 2
  array<uint8_t, UD_EDGES_POS> ud_edges_pos_pruning;
  array<uint8_t, CORNERS_POS> corners_pos_pruning;
  array<uint8_t, MIDDLE_EDGES_POS> middle_edges_pos_pruning;

  // Possible moves for phase 1 and phase 2
  const vector<Move> phase1_moves = {
      U1, U2, U3, D1, D2, D3,
      R1, R2, R3, L1, L2, L3,
      F1, F2, F3, B1, B2, B3
  };

  const vector<Move> phase2_moves = {
      U1, U2, U3, D1, D2, D3,
      R2, L2, F2, B2
  };

public:
  KociembaSolver() {
    init_arrays();
    generate_prunings_tables();
  }
private:
  void init_arrays() {
    edges_ori_pruning.fill(INIT_VAL);
    corners_ori_pruning.fill(INIT_VAL);
    middle_edges_comb_pruning.fill(INIT_VAL);

    ud_edges_pos_pruning.fill(INIT_VAL);
    corners_pos_pruning.fill(INIT_VAL);
    middle_edges_pos_pruning.fill(INIT_VAL);
  }
  void generate_prunings_tables() {
    generate_pruning_table(edges_ori_pruning, phase1_moves, &RubiksCube::get_edge_orientation_hash);
    generate_pruning_table(corners_ori_pruning, phase1_moves, &RubiksCube::get_corner_orientation_hash);
    generate_pruning_table(middle_edges_comb_pruning, phase1_moves, &RubiksCube::get_middle_comb_hash);

    generate_pruning_table(ud_edges_pos_pruning, phase2_moves, &RubiksCube::get_ud_edges_position_hash);
    generate_pruning_table(corners_pos_pruning, phase2_moves, &RubiksCube::get_corner_position_hash);
    generate_pruning_table(middle_edges_pos_pruning, phase2_moves, &RubiksCube::get_middle_edges_position_hash);
  }
  // Templated generation of all pruning tables using bfs
  template<size_t N>
  void generate_pruning_table(array<uint8_t, N>& pruning, const vector<Move>& moves, uint32_t(RubiksCube::* hash_func)() const) {
    queue<pair<RubiksCube, uint8_t>> q;
    RubiksCube solved;
    q.push({ solved, 0 });
    pruning[(solved.*hash_func)()] = 0;

    while (!q.empty()) {
      auto [cube, depth] = q.front();
      q.pop();

      for (Move move : moves) {
        RubiksCube next = cube;
        next.apply_move(move);
        uint32_t hash = (next.*hash_func)();

        if (pruning[hash] == INIT_VAL) {
          pruning[hash] = depth + 1;
          q.push({ next, depth + 1 });
        }
      }
    }
  }


  // iterative deepening search for phase 1
  // use heuristic function of minumum required quantity of moves to phase 1 state 
  bool phase1_search(RubiksCube& cube, int depth, vector<Move>& solution, Move last_move) {
    if (depth == 0) {
      return cube.in_phase1();
    }
    if (cube.in_phase1()) {
      return true;
    }

    uint32_t edge_orientation_hash = cube.get_edge_orientation_hash();
    uint32_t corner_orientation_hash = cube.get_corner_orientation_hash();
    uint32_t middle_comb_hash = cube.get_middle_comb_hash();
    uint32_t max_val = max({ edges_ori_pruning[edge_orientation_hash], corners_ori_pruning[corner_orientation_hash], 
      middle_edges_comb_pruning[middle_comb_hash] });
    if (max_val > depth) {
      return false;
    }

    for (Move move : phase1_moves) {
      if (is_redundant(last_move, move)) {
        continue;
      }
      solution.push_back(move);
      cube.apply_move(move);

      if (phase1_search(cube, depth - 1, solution, move)) {
        return true;
      }

      cube.apply_move(inverse_move(move));
      solution.pop_back();
    }

    return false;
  }

  // iterative deepening search for phase 1
    // use heuristic function of minumum required quantity of moves to solved state 
  bool phase2_search(RubiksCube& cube, int depth, vector<Move>& solution, Move last_move) {
    if (depth == 0) {
      return cube.is_solved();
    }
    if (cube.is_solved()) {
      return true;
    }

    uint32_t corners_pos_hash = cube.get_corner_position_hash();
    uint32_t ud_edges_pos_hash = cube.get_ud_edges_position_hash();
    uint32_t middle_edges_pos_hash = cube.get_middle_edges_position_hash();

    uint32_t max_val = max({
        corners_pos_pruning[corners_pos_hash],
        ud_edges_pos_pruning[ud_edges_pos_hash],
        middle_edges_pos_pruning[middle_edges_pos_hash]
      });

    if (max_val > depth) {
      return false;
    }

    for (Move move : phase2_moves) {
      if (is_redundant(last_move, move)) {
        continue;
      }

      solution.push_back(move);
      cube.apply_move(move);

      if (phase2_search(cube, depth - 1, solution, move)) {
        return true;
      }

      cube.apply_move(inverse_move(move));
      solution.pop_back();
    }

    return false;
  }

  bool is_redundant(Move prev, Move current) {
    if (prev == NO_MOVE) return false;

    // Group moves by edges (U/D/R/L/F/B)
    int prev_face = static_cast<int>(prev) / 3;
    int curr_face = static_cast<int>(current) / 3;

    // Skipping rules:
    // 1. Same layer (U after U)
    // 2. Opposite layers (U after D, R after L, etc.)
    return (prev_face == curr_face) ||
      ((prev_face / 2 == curr_face / 2) && (prev_face != curr_face));
  }

  Move inverse_move(Move move) {
    switch (move) {
    case U1: return U3;
    case U3: return U1;
    case D1: return D3;
    case D3: return D1;
    case R1: return R3;
    case R3: return R1;
    case L1: return L3;
    case L3: return L1;
    case F1: return F3;
    case F3: return F1;
    case B1: return B3;
    case B3: return B1;
    default: return move; // U2, D2, R2, L2, F2, B2
    }
  }

public:
  // Main solving method
  vector<Move> solve(RubiksCube cube) {
    vector<Move> solution;

    // Phase 1
    // 12 - maximum length of phase 1
    for (int depth = 0; depth <= 12; depth++) {
      vector<Move> phase1_solution;
      RubiksCube cube_copy = cube;

      if (phase1_search(cube_copy, depth, phase1_solution, NO_MOVE)) {
        // Applying phase 1 moves
        for (Move move : phase1_solution) {
          cube.apply_move(move);
          solution.push_back(move);
        }
        break;
      }
    }

    // Checking that cube in phase 1 state
    if (!cube.in_phase1()) {
      cerr << "Error: Phase1 did not reach target state!" << endl;
      return {};
    }

    // Phase 2
    // 18 - maximum length of phase 2
    for (int depth = 0; depth <= 18; depth++) {
      vector<Move> phase2_solution;
      RubiksCube cube_copy = cube;

      if (phase2_search(cube_copy, depth, phase2_solution, NO_MOVE)) {
        // Applying phase 2 moves
        for (Move move : phase2_solution) {
          solution.push_back(move);
        }
        break;
      }
    }

    optimize_solution(solution);

    return solution;
  }

private:
  // Optimization of the sequence of moves
  void optimize_solution(vector<Move>& solution) {
    if (solution.size() < 2) return;

    vector<Move> optimized;
    optimized.push_back(solution[0]);

    for (size_t i = 1; i < solution.size(); i++) {
      Move last = optimized.back();
      Move current = solution[i];

      // Checking if moves can be combined
      if (same_face(last, current)) {
        optimized.pop_back();
        Move combined = combine_moves(last, current);
        if (combined != NO_MOVE) {
          optimized.push_back(combined);
        }
      }
      else {
        optimized.push_back(current);
      }
    }

    solution = optimized;
  }

  bool same_face(Move a, Move b) {
    return (static_cast<int>(a) / 3) == (static_cast<int>(b) / 3);
  }

  Move combine_moves(Move a, Move b) {
    int a_base = static_cast<int>(a) % 3 + 1;
    int b_base = static_cast<int>(b) % 3 + 1;
    int face = static_cast<int>(a) / 3;

    int combined = (a_base + b_base) % 4;
    if (combined == 0) return NO_MOVE;
    return static_cast<Move>(face * 3 + combined - 1);
  }
};

// =================
// Helper functions
// =================

void print_solution(const vector<Move>& solution) {
  cout << "Solution (" << solution.size() << " moves): ";
  for (Move move : solution) {
    cout << move_names[move] << " ";
  }
  cout << endl;
}

// ==========
// CubeParser 
// ==========
// From the sequence of facelet colors in order U F R B L D from top to bottom, left to right
// To edge/corner position/orientation
class CubeParser {
private:
  static const int NUMBER_EDGE_FACELET = 2;
  static const int NUMBER_EDGES = 12;
  static const int NUMBER_CORNER_FACELET = 3;
  static const int NUMBER_CORNERS = 8;
  array<array<char, NUMBER_EDGE_FACELET>, NUMBER_EDGES> edges;
  array<array<char, NUMBER_CORNER_FACELET>, NUMBER_CORNERS> corners;

  enum EdgePositions {
    UR, UF, UL, UB, DR, DF, DL, DB, FR, FL, BL, BR
  };
  enum CornerPositions {
    URF, UFL, ULB, UBR, DFR, DLF, DBL, DRB
  };

  // Arranged the colors in a special order, the first edge color was always
  // in the priority position 1) U/D 2) L/R, for corners the first color is on U/D, the second
  // color and the third color from left to right for the top and vice versa for the bottom
  void FillEdgesCorneres(const string& state) {
    edges[UR] = { state[5], state[19] };
    edges[UF] = { state[7], state[10] };
    edges[UL] = { state[3], state[37] };
    edges[UB] = { state[1], state[28] };
    edges[DR] = { state[50], state[25] };
    edges[DL] = { state[48], state[43] };
    edges[DF] = { state[46], state[16] };
    edges[DB] = { state[52], state[34] };
    edges[FR] = { state[21], state[14] };
    edges[FL] = { state[41], state[12] };
    edges[BL] = { state[39], state[32] };
    edges[BR] = { state[23], state[30] };

    corners[URF] = { state[8], state[11], state[18] };
    corners[UFL] = { state[6], state[38], state[9] };
    corners[ULB] = { state[0], state[29], state[36] };
    corners[UBR] = { state[2], state[20], state[27] };
    corners[DFR] = { state[47], state[24], state[17] };
    corners[DLF] = { state[45], state[15], state[44] };
    corners[DBL] = { state[51], state[42], state[35] };
    corners[DRB] = { state[53], state[33], state[26] };
  }

  // Function to create a map for edge positions
  map<std::set<char>, int> createEdgePositionMap() {
    map<std::set<char>, int> map = {
        {{'U', 'R'}, 0},
        {{'U', 'F'}, 1},
        {{'U', 'L'}, 2},
        {{'U', 'B'}, 3},
        {{'D', 'R'}, 4},
        {{'D', 'F'}, 5},
        {{'D', 'L'}, 6},
        {{'D', 'B'}, 7},
        {{'F', 'R'}, 8},
        {{'F', 'L'}, 9},
        {{'B', 'L'}, 10},
        {{'B', 'R'}, 11}
    };
    return map;
  }

  // Function to create map for corner positions
  map<std::set<char>, int> createCornerPositionMap() {
    map<std::set<char>, int> map = {
        {{'U', 'R', 'F'}, 0},
        {{'U', 'F', 'L'}, 1},
        {{'U', 'L', 'B'}, 2},
        {{'U', 'B', 'R'}, 3},
        {{'D', 'F', 'R'}, 4},
        {{'D', 'L', 'F'}, 5},
        {{'D', 'B', 'L'}, 6},
        {{'D', 'R', 'B'}, 7}
    };
    return map;
  }

  bool IsUD(char c) {
    return c == 'U' || c == 'D';
  }

  bool IsLR(char c) {
    return c == 'L' || c == 'R';
  }
  
public:
  void Parse(RubiksCube& cube, const string& state) { 
    if (state.length() != 54) {
      throw invalid_argument("State string must have exactly 54 characters");
    }

    FillEdgesCorneres(state);
    auto corner_position = createCornerPositionMap();
    auto edge_position = createEdgePositionMap();
    int idx = 0;
    // If UD color on FB edges, then flip
    // Otherwise, if the color on UD or LR edges is UD LR, then they won't flip, otherwise they will flip
    for (auto& [c1, c2] : edges) {
      uint8_t flipped = 0;
      if (IsUD(c2)) {
        flipped = 1;
      }
      else if (IsUD(c1) || IsLR(c1)) {
        flipped = 0;
      }
      else {
        flipped = 1;
      }
      cube.edge_orientations[idx] = flipped;
      cube.edges[idx++] = edge_position[{c1, c2}];
    }

    idx = 0;
    // colors are arranged so that c1 is UD, c2 is left for top/right for bottom, c3 is the remaining
    for (auto& [c1, c2, c3] : corners) {
      uint8_t twist = 0;
      if (IsUD(c2)) {
        twist = 2;
      }
      else if (IsUD(c3)) {
        twist = 1;
      }

      cube.corner_orientations[idx] = twist;
      cube.corners[idx++] = corner_position[{c1, c2, c3}];
    }
  }
};

class Solver {
public:
  Solver() = default;
  void Solve(const char* state, char* solution) {
    CubeParser parser;
    parser.Parse(cube, string(state));
    
    KociembaSolver solver;
    vector<Move> moves = solver.solve(cube);

    ConvertSolution(moves, solution);
  }
private:
  RubiksCube cube;

  void ConvertSolution(const vector<Move>& moves, char* solution) {
    string str_move;
    string res;
    for (const auto& move : moves) {
      str_move = move_names[move];
      if (str_move.size() == 2 && str_move[1] == '2') {
        res += str_move[0] + string(" 0 ") + str_move[0] + string(" 0");
      }
      else if (str_move.size() == 2) {
        res += str_move[0] + string(" 0");
      }
      else {
        res += str_move[0] + string(" 1");
      }
      res += " ";
    }

    strncpy(solution, res.c_str(), res.size() + 1);
  }
};

extern "C" {
  // Creates a Solver object on the heap and returns a pointer
  Solver* Solver_new() {
    return new Solver();
  }

  // Frees memory
  void Solver_delete(Solver* solver) {
    delete solver;
  }

  // Solves a Rubik's cube and writes the solution to the buffer
  void Solver_solve(Solver* solver, const char* state, char* solution_buffer, int buffer_size) {

    // Call the Solve method, passing a temporary buffer
    solver->Solve(state, solution_buffer); // solution_buffer will be modified internally

    // If the solution does not fit into the buffer, truncate it
    if (strlen(solution_buffer) >= buffer_size - 1) {
      solution_buffer[buffer_size - 1] = '\0'; // Guarantee a null-terminated string
    }
  }
}
