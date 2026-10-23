"""Reproduce D8 / K11.4 Fig. B1 pixel comparison (no aeroelast imports).

Run from the repository root: scripts/aeroenv.sh python tools/d8_bernardi_figb1.py
Downloads are data only; hashes are verified before parsing.
"""
from pathlib import Path
import hashlib
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'logs/check/d8'
INPUTS = OUT / 'inputs'
BASE = 'https://raw.githubusercontent.com/IEAWindSystems/IEA-15-240-RWT/v1.0/OpenFAST/IEA-15-240-RWT/'
SOURCES = {
    'f14.png': ('https://wes.copernicus.org/articles/11/2345/2026/wes-11-2345-2026-f14.png', '148012b8'),
    'v1.0-blade.dat': (BASE + 'IEA-15-240-RWT_BeamDyn_blade.dat', '625253b5'),
    'v1.0-BeamDyn.dat': (BASE + 'IEA-15-240-RWT_BeamDyn.dat', 'a854b396'),
}


def download():
    INPUTS.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, (url, prefix) in SOURCES.items():
        path = INPUTS / name
        if not path.exists():
            data = urllib.request.urlopen(url, timeout=120).read()
            if not hashlib.sha256(data).hexdigest().startswith(prefix):
                raise ValueError(f'Unexpected sha256 for {url}')
            path.write_bytes(data)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if not digest.startswith(prefix):
            raise ValueError(f'Unexpected sha256 for {path}')
        hashes[name] = digest
    return hashes


def separable(v0, v1, delta):
    return bool(abs(v1 - v0) / v1 > 2 * delta)


def consensus(raw, principal):
    return raw if raw == principal else 'cannot tell'


def decide(rows):
    import numpy as np
    selected = [r for r in rows if r['separable']]
    n = len(selected)
    f1 = sum(r['closer'] == 'v1.1' for r in selected) / n if n else 0.0
    f0 = sum(r['closer'] == 'v1.0' for r in selected) / n if n else 0.0
    median = float(np.median([r['error'] / max(2 * r['delta'], .05)
                              for r in selected])) if n else float('nan')
    verdict = 'cannot tell'
    if n >= 5:
        if f1 >= .8 and median <= 1:
            verdict = 'v1.1'
        elif f0 >= .8:
            verdict = 'v1.0'
    return dict(verdict=verdict, n=n, f1=f1, f0=f0, median=median)


# Pixel calibration visually verified on the original 998x804 PNG. Coordinates
# are pixel centres, origin upper left; fit all labelled y ticks, not plot bounds.
# Panel bounds: a (87,1)-(434,268), b (637,1)-(984,268),
# c (87,420)-(434,685), d (637,420)-(984,685). All axes linear.
PANELS = {
    'a': (87, 434, [246, 171, 96, 21], [0, 50, 100, 150]),
    'b': (637, 984, [246, 201, 156, 111, 66, 21], [0, 2000, 4000, 6000, 8000, 10000]),
    'c': (87, 434, [664, 592, 520, 448], [0, 500, 1000, 1500]),
    'd': (637, 984, [636, 574, 512, 449], [0, 5, 10, 15]),
}
CURVES = {
    'EI edge': ('a', 'black', 'ei_edge', 1e9),
    'EI flap': ('a', 'teal', 'ei_flap', 1e9),
    'GJ': ('a', 'red', 'gj', 1e9),
    'inertia edge': ('b', 'black', 'm44', 1),
    'inertia flap': ('b', 'teal', 'm55', 1),
    'rhoA': ('c', 'black', 'mass_dens', 1),
    'twist': ('d', 'black', 'twist', 1),
}


def calibrate():
    import numpy as np
    result = {}
    for panel, (left, right, pixels, values) in PANELS.items():
        slope, intercept = np.polyfit(pixels, values, 1)
        residual = np.max(np.abs((np.asarray(values) - intercept) / slope - pixels))
        result[panel] = (slope, intercept, residual)
    return result


def extract(image, calibration):
    """Per-column pigment extraction, with no deck-guided selection or fitting.

    Teal and red use channel contrast (including antialiasing); black uses neutral
    dark pixels. Ring markers are centred by the median of pigment rows. Solid
    edge curves use the uppermost pigment centre (edge is above flap). Dashed
    gaps interpolate only between observed columns. Border columns use nearest
    observed pigment, with endpoint uncertainty explicitly reported.
    """
    import numpy as np
    rgb = np.asarray(image.convert('RGB')).astype(float)
    r, g, b = rgb.transpose(2, 0, 1)
    masks = {
        'red': (r - g > 35) & (r - b > 25),
        'teal': (g - r > 12) & (b - r > 4) & (r < 190),
        'black': (np.max(rgb, axis=2) - np.min(rgb, axis=2) < 18)
                 & (np.max(rgb, axis=2) < 150),
    }
    result = {}
    for name, (panel, colour, _, _) in CURVES.items():
        left, right, _, _ = PANELS[panel]
        slope, intercept, _ = calibration[panel]
        # Apply the same panel legend exclusion to every pigment, including
        # coloured fallback candidates used when the black edge overlaps them.
        panel_masks = {pigment: pixels.copy() for pigment, pixels in masks.items()}
        for pixels in panel_masks.values():
            if panel == 'a':
                pixels[:85, 305:] = False  # legend including its neutral frame
            elif panel == 'b':
                pixels[:63, 855:] = False
        mask = panel_masks[colour]
        # Short ticks project inward from top/right frame. Exclude their 5 px
        # footprints, independently of any deck values or decision outcome.
        top, bottom = (8, 264) if panel in 'ab' else (427, 681)
        xs, ys = [], []
        for x in range(left + 2, right - 5):
            pixels = np.flatnonzero(mask[top:bottom, x]) + top
            candidates = []
            if len(pixels):
                if colour == 'black':
                    pixels = pixels[pixels <= pixels[0] + 3]
                candidates.append(float(np.median(pixels)))
            if colour == 'black' and panel in 'ab':
                for pigment in (('teal', 'red') if panel == 'a' else ('teal',)):
                    cp = np.flatnonzero(panel_masks[pigment][top:bottom, x]) + top
                    if len(cp):
                        candidates.append(float(np.median(cp)))
            if not candidates:
                continue
            xs.append(x)
            ys.append(min(candidates))
        if not xs:
            raise ValueError(f'No pigment for {name}')
        x = np.arange(left, right + 1)
        y = np.interp(x, xs, ys)
        result[name] = dict(eta=(x-left)/(right-left), value=slope*y+intercept,
                            x=x, y=y, observed=len(xs),
                            gap=int(max(np.diff(xs))) if len(xs) > 1 else 0,
                            resolution=2*abs(slope))
    return result


def load_version(blade_path, axis_path):
    import numpy as np
    # The shared parser is authoritative for K, M11 and principal EI.
    from beam_reference import load_beamdyn_blade
    blade = load_beamdyn_blade(blade_path, blade_length=1.0)
    lines = blade_path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if 'DISTRIBUTED PROPERTIES' in line) + 1
    inertia = []
    i = start
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) != 1:
            i += 1
            continue
        try:
            eta = float(parts[0])
        except ValueError:
            i += 1
            continue
        j = i + 1
        while not lines[j].strip():
            j += 1
        j += 6
        while not lines[j].strip():
            j += 1
        m = np.array([[float(t) for t in lines[j+k].split()] for k in range(6)])
        inertia.append((eta, m[3, 3], m[4, 4]))
        i = j + 6
    inertia = np.asarray(inertia)
    np.testing.assert_allclose(inertia[:, 0], blade.z)
    axis_lines = axis_path.read_text().splitlines()
    j = next(i for i, line in enumerate(axis_lines) if 'initial_twist' in line) + 2
    points = []
    while j < len(axis_lines):
        tokens = axis_lines[j].split()
        if len(tokens) != 4:
            break
        try:
            row = [float(t) for t in tokens]
        except ValueError:
            break
        points.append(row)
        j += 1
    points = np.asarray(points)
    if len(points) < 3:
        raise ValueError('Missing key points')
    # Figure x is radial fraction, not swept-axis arc length.
    twist_eta = points[:, 2] / points[-1, 2]
    return dict(eta=blade.z, ei_edge=blade.ei_edge, ei_flap=blade.ei_flap,
                ei_edge_ec=blade.ei_edge_ec, ei_flap_ec=blade.ei_flap_ec,
                gj=blade.gj, mass_dens=blade.mass_dens,
                m44=inertia[:, 1], m55=inertia[:, 2],
                twist=points[:, 3], twist_eta=twist_eta)


def compare(versions, digitised, principal=False):
    import numpy as np
    rows = []
    for name, (_, _, field, scale) in CURVES.items():
        key = field + '_ec' if principal and field in ('ei_edge', 'ei_flap') else field
        eta_key = 'twist_eta' if field == 'twist' else 'eta'
        # Evaluate every native station of each version, retaining its identity;
        # identical eta values in both decks are intentionally evaluated twice.
        for version, native in versions.items():
            for eta in native[eta_key]:
                if not 0 <= eta <= 1:
                    continue
                vals = [float(np.interp(eta, deck[eta_key], deck[key]) / scale)
                        for deck in versions.values()]
                d = digitised[name]
                value = float(np.interp(eta, d['eta'], d['value']))
                delta = d['resolution'] / value if value != 0 else float('inf')
                v0, v1 = vals
                near = 'tie' if abs(value-v0) == abs(value-v1) else (
                    'v1.0' if abs(value-v0) < abs(value-v1) else 'v1.1')
                rows.append(dict(curve=name, native=version, eta=eta, value=value,
                                 v0=v0, v1=v1, delta=delta,
                                 separable=separable(v0, v1, delta), closer=near,
                                 error=abs(value-v1)/v1))
    return rows


def main():
    import numpy as np
    from PIL import Image
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    hashes = download()
    image = Image.open(INPUTS / 'f14.png')
    if image.size != (998, 804):
        raise ValueError('Unexpected image dimensions')
    v1_blade = ROOT / 'tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat'
    v1_axis = ROOT / ('tests/IEA15MW/validation papers/github-IEA-15-240-RWT/'
                      'OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat')
    for path in (v1_blade, v1_axis):
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    versions = {'v1.0': load_version(INPUTS/'v1.0-blade.dat', INPUTS/'v1.0-BeamDyn.dat'),
                'v1.1': load_version(v1_blade, v1_axis)}
    calibration = calibrate()
    digitised = extract(image, calibration)
    raw = compare(versions, digitised)
    principal = compare(versions, digitised, principal=True)
    results = {'raw': decide(raw), 'principal': decide(principal)}
    verdict = consensus(results['raw']['verdict'], results['principal']['verdict'])
    report = ['# D8 Bernardi Fig. B1 — K11.4', '', f'**Verdict: {verdict}**', '',
              'Reproduce: `scripts/aeroenv.sh python tools/d8_bernardi_figb1.py`', '',
              '## Rule (fixed before extraction)',
              'All seven plotted curves; every native station of each version within [0,1]. '
              'The other deck is linearly interpolated at each native station. Coincident '
              'native stations are retained once per version, not deduplicated.',
              'Separable: |v1.1-v1.0|/v1.1 > 2 delta_rel; delta_rel = 2 pixel data units / digitised value. '
              'Ties count for neither version. At least 5 separable stations; >=80% closer '
              'to v1.1 plus median(e/max(2 delta_rel,0.05)) <=1 gives v1.1; >=80% closer '
              'to v1.0 gives v1.0; otherwise cannot tell. Raw/principal disagreement gives cannot tell.',
              'Signed denominators are retained exactly as written in K11.4, including negative '
              'twist values (not silently replaced by absolute values). Zero digitised values '
              'have infinite relative resolution. Twist x uses key-point z/117 m, the plotted radial fraction.',
              '', '## Provenance / SHA256']
    report += [f'- `{name}`: `{digest}`' for name, digest in hashes.items()]
    report += ['', '## Calibration',
               'Visually verified linear axes. Pixel origin top left. '
               'Calibrations are fixed constants, independent of deck agreement. '
               'Residuals are tick-fit residuals, not estimates of total extraction error.',
               '|panel|x pixels (0,1)|y pixels|y values|maximum residual px|2 px data units|',
               '|---|---|---|---|---:|---:|']
    for panel, (left, right, pixels, values) in PANELS.items():
        slope, _, residual = calibration[panel]
        report.append(f'|{panel}|{left}, {right}|{pixels}|{values}|{residual:.6g}|{2*abs(slope):.6g}|')
    report += ['', '## Results', '|basis|verdict|separable|closer v1.1|closer v1.0|normalised median (<=1)|',
               '|---|---|---:|---:|---:|---:|']
    for basis, result in results.items():
        report.append(f"|{basis}|{result['verdict']}|{result['n']}|{result['f1']:.6%}|"
                      f"{result['f0']:.6%}|{result['median']:.8g}|")
    report += ['', '## Digitisation inspection',
               'Visually inspected the regenerated overlay after the shared legend-mask fix: '
               'checked both legend regions, previously stray panel-a points, axis/tick regions, '
               'and the digitised traces against all seven original curves across panels a–d. '
               'The overlay is clean of visible legend/tick contamination; extracted traces follow '
               'the figure. Deck disagreement, notably rhoA, is not an extraction adjustment target. '
               'Remaining doubts are unresolved near-zero outboard overlaps, marker-centre and '
               'antialiasing bias, interpolated dash gaps, and the nearest-pigment endpoint reads. '
               'This visual check does not prove subpixel accuracy. Panel c’s printed kg/m³ '
               'label remains inconsistent with the mandated M11 kg/m interpretation.',
               '', '## Digitisation limitations',
               '- Pigment masks use no deck values. Edge is the uppermost pigment run, including '
               'coloured pigment centres at overlaps; teal and red marker centres use median pigment rows. '
               'Missing dash columns interpolate linearly. Border columns use nearest pigment '
               '(2 px root / 6 px tip footprint); endpoints carry additional horizontal uncertainty.',
               '- Outboard stiffness/inertia curves collapse near zero: 2 pixels can exceed '
               'the plotted value, so relative uncertainty is large; no stations are removed on that account.',
               '- Panel c prints kg/m³, but comparison uses M11 in kg/m as expressly required by the addendum.',
               '- Marker radius, antialiasing and interpolated dash gaps can bias a single column. '
               'Overlay is the required visual check; tick residuals do not establish subpixel curve accuracy.',
               '', '|curve|observed columns|largest observed-column gap px|', '|---|---:|---:|']
    for name, d in digitised.items():
        report.append(f"|{name}|{d['observed']}|{d['gap']}|")
    for basis, rows in (('raw', raw), ('principal', principal)):
        report += ['', f'## {basis}: per-curve station tables',
                   'Stiffness values in GN m²; inertia in kg m; rhoA in kg/m; twist in degrees.']
        for name in CURVES:
            subset = [r for r in rows if r['curve'] == name]
            n = sum(r['separable'] for r in subset)
            report += ['', f'### {name} ({n} separable)',
                       '|native deck|eta|digitised|v1.0|v1.1|delta_rel|separable|closer-to|e/max(2 delta,5%)|',
                       '|---|---:|---:|---:|---:|---:|---|---|---:|']
            for r in subset:
                report.append(f"|{r['native']}|{r['eta']:.8g}|{r['value']:.8g}|{r['v0']:.8g}|"
                              f"{r['v1']:.8g}|{r['delta']:.8g}|{r['separable']}|{r['closer']}|"
                              f"{r['error']/max(2*r['delta'],.05):.8g}|")
    (OUT/'report.md').write_text('\n'.join(report)+'\n')
    fig, ax = plt.subplots(figsize=(12.475, 10.05), dpi=100)
    ax.imshow(image)
    for name, (panel, _, field, scale) in CURVES.items():
        left, right, _, _ = PANELS[panel]
        slope, intercept, _ = calibration[panel]
        d = digitised[name]
        from matplotlib.patches import Rectangle
        ytop, ybottom = (1, 268) if panel in 'ab' else (420, 685)
        clip = Rectangle((left, ytop), right-left, ybottom-ytop, transform=ax.transData)
        ax.plot(d['x'][::5], d['y'][::5], '.', color='#e000d0', markersize=2,
                label='digitised' if name == 'EI edge' else None, clip_path=clip)
        for version, deck in versions.items():
            eta_key = 'twist_eta' if field == 'twist' else 'eta'
            ax.plot(left+(right-left)*deck[eta_key], (deck[field]/scale-intercept)/slope,
                    color='#00a5ff' if version == 'v1.0' else '#ff9100',
                    linewidth=.9, alpha=.85,
                    label=version+' raw' if name == 'EI edge' else None)
            ax.lines[-1].set_clip_path(clip)
            if field in ('ei_edge', 'ei_flap'):
                ax.plot(left+(right-left)*deck['eta'], (deck[field+'_ec']/scale-intercept)/slope,
                        color='#00a5ff' if version == 'v1.0' else '#ff9100',
                        linewidth=.7, linestyle=':', alpha=.85)
                ax.lines[-1].set_clip_path(clip)
    ax.set_xlim(0, 998)
    ax.set_ylim(804, 0)
    ax.axis('off')
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1), ncol=3, fontsize=9)
    fig.savefig(OUT/'overlay.png', bbox_inches='tight')
    plt.close(fig)
    print(f'Verdict: {verdict}; {results}')
    print('Calibration max residual px:', {p: c[2] for p, c in calibration.items()})


def self_test():
    """Deterministic decision-rule boundaries, independent of figure/decks."""
    import numpy as np
    def rows(n, winners, error=0.01):
        return [dict(separable=True, closer='v1.1' if i < winners else 'v1.0',
                     error=error, delta=0.01) for i in range(n)]
    assert decide(rows(4, 4))['verdict'] == 'cannot tell'
    assert decide(rows(5, 4))['verdict'] == 'v1.1'
    assert decide(rows(5, 1))['verdict'] == 'v1.0'
    assert decide(rows(6, 4))['verdict'] == 'cannot tell'
    assert decide(rows(5, 5, 0.05))['verdict'] == 'v1.1'
    assert decide(rows(5, 5, 0.05001))['verdict'] == 'cannot tell'
    r = rows(5, 5, 0.08)
    for p in r:
        p['delta'] = 0.05
    assert np.isclose(decide(r)['median'], 0.8)
    assert consensus('v1.1', 'v1.0') == 'cannot tell'
    assert consensus('v1.1', 'v1.1') == 'v1.1'
    assert separable(110, 100, 0.05) is False  # strict >, not >=
    assert separable(111, 100, 0.05) is True
    # Regression: coloured legend pixels must never enter edge fallbacks.
    from PIL import Image, ImageDraw
    clean = Image.new('RGB', (998, 804), 'white')
    draw = ImageDraw.Draw(clean)
    for left, right in ((89, 428), (639, 978)):
        draw.line((left, 200, right, 200), fill=(50, 50, 50))
        draw.line((left, 210, right, 210), fill=(55, 90, 95))
    draw.line((89, 220, 428, 220), fill=(200, 20, 20))
    draw.line((89, 500, 428, 500), fill=(50, 50, 50))
    draw.line((639, 500, 978, 500), fill=(50, 50, 50))
    contaminated = clean.copy()
    draw = ImageDraw.Draw(contaminated)
    draw.line((310, 40, 420, 40), fill=(55, 90, 95))
    draw.line((310, 60, 420, 60), fill=(200, 20, 20))
    draw.line((860, 40, 970, 40), fill=(55, 90, 95))
    expected = extract(clean, calibrate())
    actual = extract(contaminated, calibrate())
    for curve in CURVES:
        np.testing.assert_array_equal(actual[curve]['value'], expected[curve]['value'],
                                      err_msg=f'Legend changed {curve}')
    print('Decision-rule and legend-mask self-tests passed')


if __name__ == '__main__':
    import sys
    if '--self-test' in sys.argv:
        self_test()
    else:
        main()
