import numpy as np, sys
sys.path.insert(0, "tests")
import test_iea15mw_s3_newmark_transient as t
case = t.transient_case.__wrapped__() if hasattr(t.transient_case, "__wrapped__") else None
print("DT", t.DT, "N", t.N_STEPS, "F1_EXACT", t.F1_EXACT)
if case:
    f=case["f1_num"]; w=2*np.pi*f
    print("f1_num", f, "spatial err %", 100*(f/t.F1_EXACT-1))
    print("newmark pred f", f*(1-(w*t.DT)**2/24), "temporal %", -100*(w*t.DT)**2/24)
