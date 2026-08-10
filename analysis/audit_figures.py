# -*- coding: utf-8 -*-
"""Every image the paper includes: is the Overleaf copy the newest build?"""
import io, os, re, datetime

TEX = 'C:/Users/stirpeg2/AppData/Local/Temp/overleaf-68d39c9d6e301aadbb376c0e/main_fix.tex'
FIG = 'C:/Users/stirpeg2/AppData/Local/Temp/overleaf-68d39c9d6e301aadbb376c0e/figures'
RES = 'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results'

s = io.open(TEX, encoding='utf8').read()
used = re.findall(r'includegraphics\[[^\]]*\]\{figures/([^}]*)\}', s)

def when(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)) if os.path.exists(p) else None

print('%-34s %-17s %-17s %s' % ('figure', 'in Overleaf', 'in results/', 'verdict'))
stale = 0
for f in sorted(set(used)):
    o = when(os.path.join(FIG, f))
    r = when(os.path.join(RES, f))
    os_ = o.strftime('%m-%d %H:%M') if o else 'MISSING'
    rs = r.strftime('%m-%d %H:%M') if r else '-- (not built here)'
    if o and r and r > o + datetime.timedelta(seconds=60):
        v = '<-- OVERLEAF IS STALE'
        stale += 1
    elif not o:
        v = '<-- MISSING'
        stale += 1
    elif not r:
        v = 'no local build'
    else:
        v = 'ok'
    print('%-34s %-17s %-17s %s' % (f, os_, rs, v))
print()
print('%d of %d figures need attention' % (stale, len(set(used))))
