#!/usr/bin/env python3
"""Tolerance-aware numerical diff for bench output .dat files.

Usage: numdiff.py [--rtol R] [--atol A] baseline.dat candidate.dat

Exit code 0 if every numeric token is within max(atol, rtol*max(|a|,|b|)),
1 otherwise. Header lines (starting with '#') and blanks must match exactly.
"""
import argparse
import sys


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-10)
    p.add_argument("baseline")
    p.add_argument("candidate")
    args = p.parse_args()

    with open(args.baseline) as fa, open(args.candidate) as fb:
        la = fa.readlines()
        lb = fb.readlines()

    if len(la) != len(lb):
        print(f"FAIL {args.candidate}: line count {len(la)} vs {len(lb)}")
        return 1

    max_rel = 0.0
    max_abs = 0.0
    for ln, (sa, sb) in enumerate(zip(la, lb), start=1):
        sa_stripped = sa.strip()
        sb_stripped = sb.strip()
        if not sa_stripped and not sb_stripped:
            continue
        if sa.startswith("#") or sb.startswith("#"):
            if sa != sb:
                print(f"FAIL {args.candidate}:{ln}: header mismatch")
                return 1
            continue
        ta, tb = sa.split(), sb.split()
        if len(ta) != len(tb):
            print(f"FAIL {args.candidate}:{ln}: column count {len(ta)} vs {len(tb)}")
            return 1
        for col, (a, b) in enumerate(zip(ta, tb), start=1):
            try:
                xa, xb = float(a), float(b)
            except ValueError:
                if a != b:
                    print(f"FAIL {args.candidate}:{ln}:col{col}: '{a}' vs '{b}'")
                    return 1
                continue
            d = abs(xa - xb)
            mag = max(abs(xa), abs(xb))
            if d > args.atol and d > args.rtol * mag:
                print(f"FAIL {args.candidate}:{ln}:col{col}: "
                      f"{xa!r} vs {xb!r} (abs={d:.3e}, rel={(d/mag if mag else 0):.3e})")
                return 1
            if mag > 0:
                max_rel = max(max_rel, d / mag)
            max_abs = max(max_abs, d)
    print(f"ok  {args.candidate}: max_rel={max_rel:.2e} max_abs={max_abs:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
