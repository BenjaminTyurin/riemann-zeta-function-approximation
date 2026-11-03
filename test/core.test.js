import test from 'node:test';
import assert from 'node:assert/strict';

import { zeta, zetaTerms, eulerMaclaurinZeta } from '../src/index.js';

/**
 * A small tolerance that reflects the accuracy we actually claim.  The
 * Euler-Maclaurin approximation with the default parameters is accurate to
 * well below 1e-12 for s in (1, 32], so 1e-9 is a comfortable bound that
 * still catches real regressions.
 */
const TOL = 1e-9;

/**
 * Known values of the Riemann zeta function at integer arguments, taken
 * from the closed forms (pi^2/6, pi^4/90, ...).  These are the ground truth
 * for the tests below.
 */
const KNOWN = {
  2: Math.PI ** 2 / 6,
  4: Math.PI ** 4 / 90,
  6: Math.PI ** 6 / 945,
  8: Math.PI ** 8 / 9450,
};

function close(got, want, tol = TOL) {
  return Math.abs(got - want) < tol;
}

test('zeta(2) matches pi^2/6', () => {
  assert.ok(close(zeta(2), KNOWN[2]));
});

test('zeta(4) matches pi^4/90', () => {
  assert.ok(close(zeta(4), KNOWN[4]));
});

test('zeta(6) matches pi^6/945', () => {
  assert.ok(close(zeta(6), KNOWN[6]));
});

test('zeta(8) matches pi^8/9450', () => {
  assert.ok(close(zeta(8), KNOWN[8]));
});

test('zeta is decreasing for s > 1', () => {
  const a = zeta(1.5);
  const b = zeta(2);
  const c = zeta(3);
  assert.ok(a > b);
  assert.ok(b > c);
});

test('zeta approaches 1 as s grows large', () => {
  assert.ok(close(zeta(100), 1, 1e-12));
});

test('zeta rejects s <= 1', () => {
  assert.throws(() => zeta(1), RangeError);
  assert.throws(() => zeta(0), RangeError);
  assert.throws(() => zeta(-2), RangeError);
});

test('zeta rejects non-finite s', () => {
  assert.throws(() => zeta(NaN), RangeError);
  assert.throws(() => zeta(Infinity), RangeError);
  assert.throws(() => zeta(-Infinity), RangeError);
  assert.throws(() => zeta('2'), RangeError);
});

test('zeta is accurate just above 1', () => {
  // zeta(1.01) is about 100.577; the Euler-Maclaurin tail handles the
  // near-pole region where the direct sum alone would need millions of terms.
  const got = zeta(1.01);
  assert.ok(close(got, 100.577, 1e-3));
});

test('eulerMaclaurinZeta with more terms is at least as accurate', () => {
  const few = eulerMaclaurinZeta(2, 8, 0);
  const many = eulerMaclaurinZeta(2, 8, 8);
  const exact = KNOWN[2];
  assert.ok(Math.abs(many - exact) <= Math.abs(few - exact));
});

test('eulerMaclaurinZeta validates its parameters', () => {
  assert.throws(() => eulerMaclaurinZeta(2, 0, 4), RangeError);
  assert.throws(() => eulerMaclaurinZeta(2, 8, -1), RangeError);
  assert.throws(() => eulerMaclaurinZeta(2, 8, 9), RangeError);
  assert.throws(() => eulerMaclaurinZeta(1, 8, 4), RangeError);
});

test('zetaTerms reports the branch used', () => {
  assert.equal(zetaTerms(2), 16);
  assert.equal(zetaTerms(40), 32);
  assert.throws(() => zetaTerms(1), RangeError);
});
