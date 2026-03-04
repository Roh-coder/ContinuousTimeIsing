
# ContinuousTimeIsing

A compact C++ implementation of a continuous-time worldline / FK-cluster
simulation. The single source file is `ConTimeIsing.cpp` and the produced
executable is `ConTimeFromJulia`.

**Requirements**

- CMake >= 3.10
- A C++14-capable compiler (GCC, Clang, MSVC)

**Build (recommended, using CMake)**

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -- -j
```

The binary will be written to `build/ConTimeIsing`.

**Quick run**

```bash
./build/ConTimeIsing --help
```

Example (reproduces the example from the source header):

```bash
./build/ConTimeIsing --layers 64,96,128 --ratio_min 0.6 --ratio_max 1.4 --ratio_n 17 \
	--burn_in 500 --sweeps_between 50 --n_configs 200 --seed 1234 --BC 1 \
	--ts_file timeseries.dat > summary.dat
```

Behavior notes

- Progress messages are printed to `stderr` (so they don't pollute the CSV).
- Summary CSV is printed to `stdout` — redirect it to a file to save results.
- Use `--ts_file` to enable a per-configuration time-series CSV (one line per
	configuration).

Files of interest

- `ConTimeIsing.cpp` — main simulation source.
- `CMakeLists.txt` — minimal build scaffold (added).

If you'd like, I can:
- Add a small example script to run a short test-config and produce sample output.
- Move sources into `src/` and add an install target.
- Add CI (GitHub Actions) and unit tests.

