"""Minimal exact integer Laurent arithmetic (standard library only).

Extracted, with no arithmetic change, from the earlier independently written
RMatrix_Review_20261026/revised/verify_exact.py. Historical file is unchanged.
Only one-variable use is needed here; the underlying class supports several.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from time import perf_counter
import json


@dataclass(frozen=True)
class Laurent:
    terms: tuple  # canonical sorted ((integer exponent tuple, integer coefficient),...)
    nv: int

    @classmethod
    def make(cls, terms, nv):
        return cls(tuple(sorted((tuple(e), c) for e, c in terms.items() if c)), nv)

    @classmethod
    def constant(cls, n, nv):
        return cls.make({(0,) * nv: n}, nv)

    def coerce(self, other):
        return other if isinstance(other, Laurent) else self.constant(other, self.nv)

    def __bool__(self):
        return bool(self.terms)

    def __add__(self, other):
        other = self.coerce(other)
        assert self.nv == other.nv
        d = dict(self.terms)
        for e, c in other.terms:
            d[e] = d.get(e, 0) + c
        return self.make(d, self.nv)

    __radd__ = __add__

    def __neg__(self):
        return self.make({e: -c for e, c in self.terms}, self.nv)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        assert self.nv == other.nv
        d = defaultdict(int)
        for e, c in self.terms:
            for f, b in other.terms:
                d[tuple(x + y for x, y in zip(e, f))] += c * b
        return self.make(d, self.nv)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int)
        if n < 0:
            assert len(self.terms) == 1 and abs(self.terms[0][1]) == 1
            e, c = self.terms[0]
            return self.make({tuple(x * n for x in e): c ** (-n)}, self.nv)
        out = self.constant(1, self.nv)
        for _ in range(n):
            out = out * self
        return out

    def evaluate(self, values):
        return sum(Fraction(c) * product(Fraction(x) ** n for x, n in zip(values, e))
                   for e, c in self.terms)

    def serial(self):
        return [{'powers': list(e), 'coefficient': str(c)} for e, c in self.terms]


def product(items):
    x = 1
    for y in items:
        x *= y
    return x


def variable(nv, k):
    return Laurent.make({tuple(int(i == k) for i in range(nv)): 1}, nv)


def accum(d, key, val):
    if val:
        d[key] = d.get(key, 0) + val
        if not d[key]:
            del d[key]


