#!/usr/bin/env bash

if ! [[ $1 =~ [0-9]+ ]]; then
    echo "usage: $0 [number of runs]" >&2
    exit 1
fi

set -e
TIMEFORMAT='%S %U'

time_compilation() {
    time g++ -std=c++20 -Wall -O3 -I../../include "$@" -o /dev/null
}

N="$1"
for (( i = 0; i < N; ++i )); do
    time_compilation without-units.cpp 2>&1
    time_compilation with-units.cpp 2>&1
    time_compilation -xc++ \
        <(sed '1s/^/#include <iostream>\n/' without-units.cpp) \
        2>&1
done | python -c '''
import sys, math
ts = [ sum(float(t) for t in line.split()) for line in sys.stdin ]
def mean_stdev(xs):
    n = len(xs)
    m = sum(xs)/n
    s = math.sqrt(sum( (x - m)**2 for x in xs ) / (n - 1))
    return (m, s)
n = 3
t = [ mean_stdev(ts[i::n]) for i in range(n) ]
dt = [
    (t[i][0] - t[0][0], math.sqrt(t[i][1]**2 + t[0][1]**2))
    for i in range(1,n)
]
def p(m,s):
    return f"{m:.6f} ± {s:.6f}"
print(f"""\
Without visit   : { p(*t[0]) }
With    visit   : { p(*t[1]) }
Difference      : { p(*dt[0]) }
Include iostream: { p(*t[2]) }
Difference      : { p(*dt[1]) }\
""")
'''
