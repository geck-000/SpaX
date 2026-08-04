#!/usr/bin/env python3
"""Locate the hardest-to-read prose in a LaTeX manuscript.

A reviewer asked for shorter sentences and fewer parenthetical insertions set
off by em dashes. Rewriting blind risks touching the passages that are already
fine while missing the ones that prompted the request, so this measures first:
per section, the number of sentences, their mean length, the count of em-dash
insertions, and how many sentences exceed forty words.

Figures and tables are stripped before counting, since captions are read
differently from body text and their length is not the complaint.

    python3 prose_audit.py main_rev.tex
"""
import re
import sys


def body(t):
    t = re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', '', t, flags=re.S)
    t = re.sub(r'\\begin\{table\*?\}.*?\\end\{table\*?\}', '', t, flags=re.S)
    return t


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'main_rev.tex'
    s = open(path, encoding='utf8', newline='').read()
    pat = re.compile(r'\\(?:sub)*section\{([^}]*)\}')
    secs = [(m.start(), m.group(1)) for m in pat.finditer(s)]

    print('%-50s %6s %8s %8s %7s' % ('section', 'sents', 'avg wds', 'emdash', '>40w'))
    tot = 0
    rows = []
    for i, (a, name) in enumerate(secs):
        b = secs[i + 1][0] if i + 1 < len(secs) else len(s)
        t = body(s[a:b])
        d = t.count('---')
        tot += d
        sents = [x for x in re.split(r'(?<=[.!?])\s+', re.sub(r'\s+', ' ', t))
                 if len(x.split()) > 4]
        if len(sents) < 3:
            continue
        avg = sum(len(x.split()) for x in sents) / len(sents)
        l40 = sum(1 for x in sents if len(x.split()) > 40)
        rows.append((avg, name, len(sents), d, l40))
        print('%-50s %6d %8.1f %8d %7d' % (name[:50], len(sents), avg, d, l40))

    print()
    print('total em-dash insertions: %d' % tot)
    print()
    print('worst by mean sentence length:')
    for avg, name, n, d, l40 in sorted(rows, reverse=True)[:8]:
        print('  %-44s %5.1f wds  %2d dashes  %2d sents >40 wds'
              % (name[:44], avg, d, l40))


if __name__ == '__main__':
    main()
