import sys, yaml, numpy as np
a = yaml.safe_load(open(sys.argv[1])); b = yaml.safe_load(open(sys.argv[2]))
def walk(x, y, p, out):
    if isinstance(x, dict) and isinstance(y, dict):
        for k in sorted(set(x) | set(y), key=str):
            if k not in x: out.append(f"+B {p}/{k}")
            elif k not in y: out.append(f"-A {p}/{k}")
            else: walk(x[k], y[k], f"{p}/{k}", out)
    elif isinstance(x, list) and isinstance(y, list):
        if all(isinstance(v, (int, float)) for v in x + y) and x and y:
            if len(x) != len(y): out.append(f"len {p}: {len(x)} vs {len(y)}")
            else:
                xa, ya = np.array(x, float), np.array(y, float)
                d = np.max(np.abs(xa - ya) / np.maximum(np.max(np.abs(xa)), 1e-30))
                if d > 1e-9: out.append(f"num {p}: max rel {d:.3g}")
        elif all(isinstance(v, dict) and 'name' in v for v in x + y):
            dx = {v['name']: v for v in x}; dy = {v['name']: v for v in y}
            walk(dx, dy, p, out)
        elif len(x) != len(y): out.append(f"len {p}: {len(x)} vs {len(y)}")
        else:
            for i, (u, v) in enumerate(zip(x, y)): walk(u, v, f"{p}[{i}]", out)
    elif x != y:
        if isinstance(x, (int, float)) and isinstance(y, (int, float)) and abs(x - y) <= 1e-9 * max(abs(x), 1e-30): return
        out.append(f"val {p}: {str(x)[:60]} vs {str(y)[:60]}")
out = []
walk(a['components']['blade'], b['components']['blade'], 'blade', out)
used = set()
for L in a['components']['blade']['structure']['layers']: used.add(L['material'])
ma = {m['name']: m for m in a['materials']}; mb = {m['name']: m for m in b['materials']}
mo = []
for n in sorted(used): walk(ma.get(n, {}), mb.get(n, {}), f"mat/{n}", mo)
print(len(out), "blade diffs"); print("\n".join(out[:60]))
print(len(mo), "used-material diffs"); print("\n".join(mo[:30]))
