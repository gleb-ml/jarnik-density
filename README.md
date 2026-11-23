# Central symmetry of minimum-perimeter convex lattice polygons

Gleb Golovkov, 2026

## Overview

For each even integer $n\ge 4$, the paper considers whether some convex lattice $n$-gon of minimum Euclidean perimeter lacks central symmetry. It proves that the set of such side counts has a relative natural density $\delta$ among the even integers. The density satisfies $16/455\le\delta<1$.

This repository reproduces the paper's computational tables by enumerating primitive lattice circles and their zero-sum relations.

## Requirements

Python 3.8 or later. The script uses only the Python standard library.

## Usage

Run the default computation from the repository root.

```sh
python jarnik_verify.py
```

The program enumerates squared radii through `1_000_000`, displays the first 32 minimal radii, and prints the total and the four summary rows. A smaller run can be configured as follows.

```sh
python jarnik_verify.py --limit 100000 --max-k 4 --seed-limit 10000 --show-cores 10
```

`--max-k` caps the number of distinct prime divisors congruent to $1$ modulo $4$. `--seed-limit` sets the largest minimal radius used to certify strict descendants. Set `--show-cores 0` to omit the list of minimal radii.

## Reproducing the Paper Tables

The first minimal radius is $65=5\cdot13$, with $t=8$ opposite pairs and minimum relation size $s=6$. The default run found 4,603 minimal radii and took about 12 seconds in the verification environment with Python 3.12.14.

```text
 Cutoff m  Max even n  Even n >= 4   Marked  Marked %  Strict inh.  Strict inh. %       S_2       S_3
---------  ----------  -----------  -------  --------  -----------  -------------  --------  --------
   30,000      57,304       28,651    2,535   8.8479%        1,992        6.9526%  0.360481  0.043183
  100,000     191,040       95,519   10,942  11.4553%        6,976        7.3033%  0.362278  0.122565
  300,000     572,944      286,471   37,359  13.0411%       20,268        7.0751%  0.362734  0.216522
1,000,000   1,909,776      954,887  123,653  12.9495%       68,969        7.2227%  0.363093  0.305351
```

`Max even n` is $U(X)$ in the paper, and the denominator of each percentage is $U(X)/2-1$. `Marked` counts even side numbers admitting a minimizer without central symmetry. `Strict inh.` counts those certified by a relation from a proper unitary divisor in the seed list. The `S_k` columns give the partial weight sums, rounded to six decimal places.

## License

The code and documentation are available under the [MIT License](LICENSE).
