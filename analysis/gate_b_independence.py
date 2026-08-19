import glob, sys, math
sys.path.insert(0, "/scratch/project_2019020/test_rve")
from SpaX_PostProcess import extract_first_order
E_ICE, DRAIN = 9.37, 1.04
B = {"b010": 0.10, "b020": 0.20, "b030": 0.30}
rows = []
for o in sorted(glob.glob("Job-LB_p150_b0*_drn_*-utx.odb")):
    tag = o[:-4].replace("Job-", "").replace("-utx", "")
    bkey = tag.split("_")[2]
    try:
        r = extract_first_order(o, "S11", 0.005/0.5, 0.5)
        E = float(r.get("E_eff") or 0)/1e9
        phi = r.get("phi_inclusion") or 0
    except Exception as e:
        print("%-26s EXTRACT FAIL %s" % (tag, e)); continue
    if E <= 0: 
        print("%-26s E=0 (incomplete)" % tag); continue
    b = B[bkey]
    ep = E_ICE*(1-1.65*phi)/DRAIN
    n = math.log(E/ep)/math.log(b)
    rows.append((tag, b, phi, E, n))
    print("%-26s b=%.2f phi=%.4f E=%6.3f  n=%.4f" % (tag, b, phi, E, n))
if rows:
    import collections
    d = collections.defaultdict(list)
    for t,b,p,E,n in rows: d[b].append(n)
    print("")
    print("n by b at phi~0.15 (fixed phi, b swept):")
    ns=[]
    for b in sorted(d):
        m=sum(d[b])/len(d[b]); ns.append(m)
        print("  b=%.2f  n=%.4f  (%d seeds)" % (b, m, len(d[b])))
    print("")
    print("spread across b: %.4f" % (max(ns)-min(ns)))
    print("VERDICT: %s" % ("n IS b-independent -- pooling is legitimate"
          if max(ns)-min(ns) < 0.10 else
          "n DEPENDS on b -- the phi-dependence may be b in disguise"))
