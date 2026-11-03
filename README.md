# riemann-zeta-function-approximation

Approximates the Riemann zeta function for real arguments greater than one
using the Euler-Maclaurin summation formula.

## Usage

```js
import { zeta, eulerMaclaurinZeta, zetaTerms } from './src/index.js';

zeta(2);                    // 1.6449340668482264 (pi^2/6)
zeta(4);                    // 1.0823232337111382 (pi^4/90)
zeta(1.01);                 // ~100.577

eulerMaclaurinZeta(2, 8, 4); // explicit truncation: 8 direct terms, 4 Bernoulli terms
zetaTerms(2);                // 16
zetaTerms(40);               // 32
```

## Why this exists

The defining series `sum n^-s` converges for `s > 1`, but near `s = 1` it
converges so slowly that summing it directly is impractical.  The
Euler-Maclaurin formula splits the sum into a short direct part plus an
integral and a Bernoulli-number correction, which restores accuracy with only
a handful of terms.  The trade-off is that the correction is an asymptotic
series: adding more Bernoulli terms eventually makes the result worse, not
better.  This library fixes the number of correction terms at eight, which is
enough for double precision over the whole supported domain, and exposes
`eulerMaclaurinZeta` for callers who want to vary that number themselves.

## The awkward edge

The function is only defined for `s > 1`.  At `s = 1` the zeta function has a
pole, and for `s <= 1` the Euler-Maclaurin tail used here is not valid.  Any
call with `s <= 1`, `NaN`, or `Infinity` throws a `RangeError` rather than
returning a silently wrong number.  The other edge is the asymptotic nature of
the correction: `eulerMaclaurinZeta(s, n, terms)` will happily accept up to
eight Bernoulli terms, but passing more than the default is not guaranteed to
improve accuracy for every `s`.
