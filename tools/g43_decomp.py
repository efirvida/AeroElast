import numpy as np, sys
sys.path.insert(0, "tests")
import test_iea15mw_s3_newmark_transient as t
case = t.transient_case.__wrapped__() if hasattr(t.transient_case, "__wrapped__") else None
print("DT", t.DT, "N", t.N_STEPS, "F1_EXACT", t.F1_EXACT)
if case:
    f=case["f1_num"]; w=2*np.pi*f
    print("f1_num", f, "spatial err %", 100*(f/t.F1_EXACT-1))
    fd=(2/t.DT)*np.arctan(w*t.DT/2)/(2*np.pi)  # exact trapezoidal dispersion, ~ f(1-(w dt)^2/12)
    print("newmark pred f", fd, "temporal %", 100*(fd/f-1))
