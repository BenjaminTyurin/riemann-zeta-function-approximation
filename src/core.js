/**
 * Approximation of the Riemann zeta function for real arguments s > 1.
 *
 * The series sum_{n=1}^inf n^-s converges for s > 1 but converges slowly
 * near s = 1, so we use the Euler-Maclaurin summation formula to accelerate
 * it.  The formula replaces the tail of the sum by an integral plus a
 * correction built from Bernoulli numbers.  For the range of s handled here
 * (s > 1, and in practice s >= 1.01) a small number of terms is enough to
 * reach near double-precision accuracy.
 *
 * Design decision: we accept only real s > 1.  The Euler-Maclaurin formula
 * also applies to complex arguments, but implementing the complex branch
 * would double the surface area and the tests.  A caller who needs complex
 * arguments can use a different library; this one keeps a single, well-tested
 * real path.
 */

/**
 * Bernoulli numbers B_2, B_4, ... B_2k up to k = 8.
 *
 * Only even-index Bernoulli numbers appear in the Euler-Maclaurin tail, and
 * the odd ones (except B_1) are zero.  We hard-code the first eight even
 * Bernoulli numbers because computing them at runtime would need either a
 * recurrence or a table anyway, and the constants are exact rationals with
 * small denominators.  B_2k is stored as its signed rational value in
 * floating point.
 */
const BERNOULLI_EVEN = [
  1 / 6,          // B_2
  -1 / 30,        // B_4
  1 / 42,         // B_6
  -1 / 30,        // B_8
  5 / 66,         // B_10
  -691 / 2730,    // B_12
  7 / 6,          // B_14
  -3617 / 510,    // B_16
];

/**
 * The largest s for which the Euler-Maclaurin correction is still useful.
 * Above this value the direct sum converges so fast that adding correction
 * terms only introduces rounding noise.
 */
const DIRECT_SUM_CUTOFF = 32;

/**
 * The number of leading terms summed directly before the Euler-Maclaurin
 * tail is applied.  N must be at least 1.  A value around 16 balances the
 * cost of the direct sum against the accuracy of the tail expansion.
 */
const DEFAULT_N = 16;

/**
 * The number of Bernoulli correction terms used in the tail.  Each term is
 * O(N^{-2k-1}), so eight terms already push the truncation error below
 * double precision for N = 16 and s > 1.
 */
const DEFAULT_TERMS = 8;

/**
 * Validate that s is a finite real greater than 1.
 *
 * @param {number} s - the argument of the zeta function.
 * @returns {void}
 * @throws {RangeError} if s is not a finite number greater than 1.
 */
function assertValidS(s) {
  if (typeof s !== 'number' || !Number.isFinite(s)) {
    throw new RangeError('s must be a finite number');
  }
  if (s <= 1) {
    throw new RangeError('s must be greater than 1');
  }
}

/**
 * Validate the truncation parameters.
 *
 * @param {number} n - number of direct terms.
 * @param {number} terms - number of Bernoulli correction terms.
 * @returns {void}
 * @throws {RangeError} if n < 1 or terms < 0.
 */
function assertValidParams(n, terms) {
  if (!Number.isInteger(n) || n < 1) {
    throw new RangeError('n must be a positive integer');
  }
  if (!Number.isInteger(terms) || terms < 0) {
    throw new RangeError('terms must be a non-negative integer');
  }
  if (terms > BERNOULLI_EVEN.length) {
    throw new RangeError(`terms must be at most ${BERNOULLI_EVEN.length}`);
  }
}

/**
 * Compute the direct partial sum sum_{k=1}^{n} k^-s.
 *
 * @param {number} s - the argument, s > 1.
 * @param {number} n - the number of terms, n >= 1.
 * @returns {number} the partial sum.
 */
function partialSum(s, n) {
  let sum = 0;
  for (let k = 1; k <= n; k += 1) {
    sum += Math.pow(k, -s);
  }
  return sum;
}

/**
 * Apply the Euler-Maclaurin correction to a partial sum.
 *
 * The formula is:
 *   zeta(s) = sum_{k=1}^{n-1} k^-s
 *           + n^(1-s)/(s-1)
 *           + 1/2 * n^-s
 *           + sum_{j=1}^{m} B_{2j}/(2j)! * risingFactorial(s, 2j-1) * n^(-s-2j+1)
 *           + R
 * where R is the remainder.  We use the first m even Bernoulli numbers.
 *
 * @param {number} s - the argument, s > 1.
 * @param {number} n - the number of direct terms.
 * @param {number} terms - the number of Bernoulli correction terms.
 * @returns {number} the approximation.
 */
function eulerMaclaurinTail(s, n, terms) {
  let result = 0;
  result += Math.pow(n, 1 - s) / (s - 1);
  result += 0.5 * Math.pow(n, -s);

  let rising = s;
  let factorial = 2;
  let power = Math.pow(n, -s - 1);
  for (let j = 0; j < terms; j += 1) {
    const b = BERNOULLI_EVEN[j];
    result += (b / factorial) * rising * power;

    // Advance risingFactorial(s, 2j-1) to risingFactorial(s, 2j+1):
    // multiply by (s + 2j - 1)(s + 2j).
    const a = s + 2 * j;
    rising *= a * (a + 1);

    // Advance (2j)! to (2j+2)!.
    factorial *= (2 * j + 2) * (2 * j + 3);

    // Advance n^(-s-2j+1) to n^(-s-2j-1).
    power /= n * n;
  }
  return result;
}

/**
 * Approximate the Riemann zeta function for real s > 1.
 *
 * For very large s the direct sum converges so quickly that the correction
 * terms would only add rounding noise, so we switch to a pure direct sum
 * above DIRECT_SUM_CUTOFF.  This keeps the function accurate across the
 * whole domain with a single code path.
 *
 * @param {number} s - the argument of the zeta function, s > 1.
 * @returns {number} an approximation of zeta(s).
 * @throws {RangeError} if s is not a finite number greater than 1.
 */
export function zeta(s) {
  assertValidS(s);

  if (s > DIRECT_SUM_CUTOFF) {
    // The tail beyond n = 32 is smaller than double precision for s > 32,
    // so summing 32 terms is already exact to machine precision.
    return partialSum(s, DIRECT_SUM_CUTOFF);
  }

  const n = DEFAULT_N;
  const head = partialSum(s, n - 1);
  const tail = eulerMaclaurinTail(s, n, DEFAULT_TERMS);
  return head + tail;
}

/**
 * Expose the Euler-Maclaurin approximation with explicit truncation
 * parameters, for callers who want to study the convergence of the formula.
 *
 * @param {number} s - the argument, s > 1.
 * @param {number} n - the number of direct terms, n >= 1.
 * @param {number} terms - the number of Bernoulli correction terms, 0 <= terms <= 8.
 * @returns {number} the approximation.
 * @throws {RangeError} if any argument is out of range.
 */
export function eulerMaclaurinZeta(s, n, terms) {
  assertValidS(s);
  assertValidParams(n, terms);
  const head = partialSum(s, n - 1);
  const tail = eulerMaclaurinTail(s, n, terms);
  return head + tail;
}

/**
 * Return the number of terms used internally for a given s.  This is a small
 * convenience for callers who want to know which branch `zeta` took.
 *
 * @param {number} s - the argument, s > 1.
 * @returns {number} the number of direct terms.
 * @throws {RangeError} if s is not a finite number greater than 1.
 */
export function zetaTerms(s) {
  assertValidS(s);
  return s > DIRECT_SUM_CUTOFF ? DIRECT_SUM_CUTOFF : DEFAULT_N;
}
