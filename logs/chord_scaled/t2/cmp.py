import json, sys, numpy as np
r = json.load(open(sys.argv[1])); deck = sys.argv[2]
lines = [l.split() for l in open(deck).read().split("DISTRIBUTED PROPERTIES")[1].split("\n")[1:]]
nums = [l for l in lines if l]
st = {}; i = 0
while i < len(nums):
    if len(nums[i]) == 1:
        eta = float(nums[i][0]); rows = [r for r in nums[i+1:i+7]]
        assert all(len(r) == 6 for r in rows), (eta, rows)
        st[round(eta, 4)] = np.array(rows, float); i += 13
    else: i += 1
print("deck K66 at 0.5:", st[0.5][5, 5], " diag:", np.diag(st[0.5]))
def cond(K, t, s):
    s = list(s); return K[t, t] - K[t, s] @ np.linalg.solve(K[np.ix_(s, s)], K[s, t])
print(" eta  | EA S/D | GJ S/D raw | GJ cond-shear S/D | EI_a S/D | EI_b S/D | GAa S/D | GAb S/D")
for x, TS in zip(r["stations"], r["TS"]):
    T = np.array(TS); D = st[round(x, 4)]
    gs = cond(T, 3, [1, 2]); gd = cond(D, 5, [0, 1])
    print(f" {x:.2f} | {T[0,0]/D[2,2]:.3f} | {T[3,3]/D[5,5]:.3f} | {gs/gd:.3f} ({gs:.4e} vs {gd:.4e}) | "
          f"{T[4,4]/D[3,3]:.3f} {T[4,4]/D[4,4]:.3f} | {T[5,5]/D[4,4]:.3f} {T[5,5]/D[3,3]:.3f} | {T[1,1]/D[0,0]:.3f} {T[1,1]/D[1,1]:.3f} | {T[2,2]/D[1,1]:.3f} {T[2,2]/D[0,0]:.3f}")
print("sha", r["sha256"][:16])
print("invariant bending (Schur onto bending pair after axial, eigenvalues) S vs D:")
for x, TS in zip(r["stations"], r["TS"]):
    T = np.array(TS); D = st[round(x, 4)]
    def bend(K, ax, b):
        idx = [ax] + b; S = K[np.ix_(idx, idx)]
        B = S[1:, 1:] - np.outer(S[1:, 0], S[0, 1:]) / S[0, 0]
        return np.sort(np.linalg.eigvalsh(B))
    bs = bend(T, 0, [4, 5]); bd = bend(D, 2, [3, 4])
    print(f" {x:.2f}  EI_min {bs[0]:.4e}/{bd[0]:.4e} = {bs[0]/bd[0]:.3f}   EI_max {bs[1]:.4e}/{bd[1]:.4e} = {bs[1]/bd[1]:.3f}")
