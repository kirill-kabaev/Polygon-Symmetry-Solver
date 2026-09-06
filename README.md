# 2D Convex Polygon Symmetry Detector

Technical assessment implementation for **Shenzhen Honeycomb Interconnection Technology Co., Ltd. (Geometry Technology Department)**.

- **Author:** Kabaev Kirill
---

## Features & Compliance

- **C++ Standard:** C++17 (zero external dependencies, standard library only).
- **Documentation:** Full Doxygen documentation for all structures and methods.
- **Robustness:** Precision tolerance $\varepsilon = 10^{-5}$ with cyclic orientation checks.
- **Complexity:** $\mathcal{O}(n^2)$ time, $\mathcal{O}(n)$ memory.

---

## How It Works

1. **Centroid:** Calculates polygon center of mass $C$. Any symmetry axis must pass through $C$.
2. **Candidates:** Pairs $C$ with each vertex and edge midpoint ($2n$ features), filtering duplicate angles in $[0, \pi)$.
3. **Validation:** Reflects all vertices across candidate lines and verifies metric matching and reverse cyclic ordering.
4. **Output:** Extracts the two extreme boundary points on opposite sides of $C$.

---

## Build

```bash
# Linux / macOS (GCC or Clang)
g++ -O3 -std=c++17 symmetry_test.cpp -o symmetry_detector

# Windows (MSVC)
cl /O2 /std:c++17 symmetry_test.cpp /Fe:symmetry_detector.exe
```

## Usage

```bash
# Run with input file
./symmetry_detector polygon.txt

# Run built-in sample tests
./symmetry_detector
```

