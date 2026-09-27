"""Connect LED junction nets with solid copper polygons (Draw Polygon shapes)
on F.Cu + B.Cu via the KiCad IPC API.  v4: octilinear geometry.

Shape = octilinear hull (all edges at 0/45/90/135 deg) of the junction's pad
bounding boxes (+margin), minus octilinear clearance cuts around other-net
copper and the board edge, then de-slivered with a mitre-join opening that
preserves the octilinear edge directions.

Pad geometry is taken from a pcbnew dump of the live board (exported via the
IPC API) because get_pad_shapes_as_polygons returns offset shapes for
rotated/fractionally-placed footprints (KiCad 10.0.5 IPC bug).  Pads with no
net inherit the net of a same-numbered sibling pad on the same footprint.

Usage: connect_leds_poly.py NET_NAME [NET_NAME ...]
       connect_leds_poly.py --all
       connect_leds_poly.py --redo        (delete existing junction polys, redo all)
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import time

from kipy import KiCad
import kipy.board_types as bt
from kipy.geometry import PolygonWithHoles, PolyLine, PolyLineNode
from kipy.util.units import from_mm, to_mm
from shapely.geometry import Polygon, MultiPolygon, LineString, Point, MultiPoint, box
from shapely.ops import unary_union

MARGINS = [from_mm(0.15), from_mm(0.4), from_mm(0.8), from_mm(1.2)]
CLEARANCE = from_mm(0.2)         # netclass 'Default' clearance
CLEAR_EPS = from_mm(0.025)       # safety on top of clearance
EDGE_CLEARANCE = from_mm(0.55)   # board rule min_copper_edge_clearance 0.5 + safety
SMOOTH_R = from_mm(0.2)          # sliver/spike removal radius (mitre opening)
SIMPLIFY = from_mm(0.001)        # collinear-vertex cleanup only
ARC_STEP = from_mm(0.2)          # Edge.Cuts arc discretization
LED_FP = 'IN-P36ATEUW'
F_CU = bt.BoardLayer.Value('BL_F_Cu')
B_CU = bt.BoardLayer.Value('BL_B_Cu')
HELPER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'pad_dump_helper.py')


def retry(fn, *a, **kw):
    for i in range(20):
        try:
            return fn(*a, **kw)
        except Exception as e:
            if 'busy' in str(e) and i < 19:
                time.sleep(0.5)
                continue
            raise


def pwh_to_shapely(pwh):
    shell = [(n.point.x, n.point.y) for n in pwh.outline.nodes]
    holes = [[(n.point.x, n.point.y) for n in h.nodes] for h in pwh.holes]
    return Polygon(shell, holes)


def shapely_to_pwh(poly):
    pwh = PolygonWithHoles()
    for x, y in list(poly.exterior.coords)[:-1]:
        pwh.outline.append(PolyLineNode.from_xy(int(round(x)), int(round(y))))
    pwh.outline.closed = True
    for interior in poly.interiors:
        h = PolyLine()
        for x, y in list(interior.coords)[:-1]:
            h.append(PolyLineNode.from_xy(int(round(x)), int(round(y))))
        h.closed = True
        pwh.holes.append(h)
    return pwh


def arc_points(arc):
    x1, y1 = arc.start.x, arc.start.y
    x2, y2 = arc.mid.x, arc.mid.y
    x3, y3 = arc.end.x, arc.end.y
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-9:
        return [(x1, y1), (x3, y3)]
    ux = ((x1**2 + y1**2) * (y2 - y3) + (x2**2 + y2**2) * (y3 - y1)
          + (x3**2 + y3**2) * (y1 - y2)) / d
    uy = ((x1**2 + y1**2) * (x3 - x2) + (x2**2 + y2**2) * (x1 - x3)
          + (x3**2 + y3**2) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    a1 = math.atan2(y1 - uy, x1 - ux)
    am = math.atan2(y2 - uy, x2 - ux)
    a3 = math.atan2(y3 - uy, x3 - ux)
    ccw = (a3 - a1) % (2 * math.pi)
    amid_ccw = (am - a1) % (2 * math.pi)
    sweep = ccw if amid_ccw <= ccw else ccw - 2 * math.pi
    n = max(4, int(abs(sweep) * r / ARC_STEP))
    return [(ux + r * math.cos(a1 + sweep * i / n),
             uy + r * math.sin(a1 + sweep * i / n)) for i in range(n + 1)]


def board_area(board):
    segs = []
    for s in retry(board.get_shapes):
        if s.layer != bt.BoardLayer.Value('BL_Edge_Cuts'):
            continue
        if isinstance(s, bt.BoardArc):
            segs.append(arc_points(s))
        elif isinstance(s, bt.BoardSegment):
            segs.append([(s.start.x, s.start.y), (s.end.x, s.end.y)])
        elif isinstance(s, bt.BoardRectangle):
            tl, br = s.top_left, s.bottom_right
            segs.append([(tl.x, tl.y), (br.x, tl.y), (br.x, br.y),
                         (tl.x, br.y), (tl.x, tl.y)])
        elif isinstance(s, bt.BoardCircle):
            c, r = s.center, (s.radius_point - s.center).length()
            segs.append([(c.x + r * math.cos(t), c.y + r * math.sin(t))
                         for t in [i / 64 * 2 * math.pi for i in range(65)]])
    rings = [s for s in segs if s[0] == s[-1] and len(s) > 3]
    open_segs = [s for s in segs if s[0] != s[-1]]
    TOL = from_mm(0.01)
    while open_segs:
        chain = list(open_segs.pop())
        progress = True
        while progress and math.hypot(chain[0][0] - chain[-1][0],
                                      chain[0][1] - chain[-1][1]) > TOL:
            progress = False
            for i, s in enumerate(open_segs):
                for pts in (s, s[::-1]):
                    if math.hypot(pts[0][0] - chain[-1][0],
                                  pts[0][1] - chain[-1][1]) <= TOL:
                        chain += pts[1:]
                        open_segs.pop(i)
                        progress = True
                        break
                if progress:
                    break
        rings.append(chain)
    polys = sorted((Polygon(r).buffer(0) for r in rings), key=lambda p: -p.area)
    area = polys[0]
    for p in polys[1:]:
        area = area.difference(p) if area.contains(p) else area.union(p)
    return area


def dump_pads(board):
    """Export the live board and get exact pad geometry via pcbnew."""
    with tempfile.NamedTemporaryFile('w', suffix='.kicad_pcb',
                                     delete=False) as f:
        f.write(retry(board.get_as_string))
        path = f.name
    try:
        res = subprocess.run(['python3', HELPER, path],
                             capture_output=True, text=True, check=True)
        pads = json.loads(res.stdout)
    finally:
        os.unlink(path)
    # effective net: inherit from same-numbered sibling on same footprint
    num_net = {}
    for p in pads:
        if p['net']:
            num_net.setdefault((p['ref'], p['num']), p['net'])
    for p in pads:
        p['eff'] = p['net'] or num_net.get((p['ref'], p['num']), '')
        g = Polygon(p['poly']) if len(p['poly']) >= 3 else None
        p['box'] = box(*g.bounds) if g is not None else None
    return pads


SQ2 = math.sqrt(2)


def octo_hull(points, grow=0.0):
    """Smallest octilinear (0/45/90/135 deg edges) convex polygon containing
    the points, grown outward by `grow` in every direction."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    us = [p[0] + p[1] for p in points]
    vs = [p[0] - p[1] for p in points]
    g2 = grow * SQ2
    bx = box(min(xs) - grow, min(ys) - grow, max(xs) + grow, max(ys) + grow)
    B = 1e10  # half-planes as very large boxes, rotated 45 deg for u/v
    umin, umax = min(us) - g2, max(us) + g2
    ustrip = Polygon([(umin - B, B), (umin + B, -B),
                      (umax + B, -B), (umax - B, B)])
    vmin, vmax = min(vs) - g2, max(vs) + g2
    vstrip = Polygon([(vmin - B, -B), (vmin + B, B),
                      (vmax + B, B), (vmax - B, -B)])
    return bx.intersection(ustrip).intersection(vstrip)


def obstacle_oct(geom, grow):
    """Octilinear over-approximation of an obstacle, grown by clearance."""
    geoms = geom.geoms if hasattr(geom, 'geoms') else [geom]
    parts = []
    for g in geoms:
        if g.is_empty:
            continue
        if isinstance(g, Point):
            parts.append(octo_hull([(g.x, g.y)], grow))
        elif isinstance(g, Polygon):
            parts.append(octo_hull(list(g.exterior.coords), grow))
        else:
            parts.append(octo_hull(list(g.coords), grow))
    return unary_union(parts)


def smooth(poly, end_union, refs):
    """Remove slivers/spikes with a mitre-join opening, which preserves
    octilinear edge directions.  Mitre dilation can overshoot at corners
    created by the erosion, so clip back to the original shape."""
    s = (poly.buffer(-SMOOTH_R, join_style=2)
             .buffer(SMOOTH_R, join_style=2)
             .intersection(poly))
    pieces = list(s.geoms) if isinstance(s, MultiPolygon) else [s]
    keep = [pc for pc in pieces if not pc.is_empty
            and all(pc.intersects(end_union[r]) for r in refs)]
    if keep:
        return max(keep, key=lambda pc: pc.area).simplify(SIMPLIFY)
    return poly.simplify(SIMPLIFY)


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 1

    k = KiCad()
    b = k.get_board()

    pads = dump_pads(b)
    led_refs = {p['ref'] for p in pads if LED_FP in p['fpid']}
    nets_by_name = {n.name: n for n in retry(b.get_nets)}

    junction_pads = {}
    for p in pads:
        if p['ref'] in led_refs and p['eff'] and p['box'] is not None:
            junction_pads.setdefault(p['eff'], []).append(p)

    junction_nets = {n for n, ps in junction_pads.items()
                     if len({p['ref'] for p in ps}) == 2
                     and all(p['ref'] in led_refs for p in ps)}

    redo = args == ['--redo']
    if args in (['--all'], ['--redo']):
        targets = sorted(junction_nets)
    else:
        targets = args

    # explicit net names mean "replace": drop their old polygons first
    doomed_nets = junction_nets if redo else set(targets)
    shapes = retry(b.get_shapes)
    old = [s for s in shapes if isinstance(s, bt.BoardPolygon)
           and s.layer in (F_CU, B_CU) and s.net.name in doomed_nets]
    if old:
        print(f'Deleting {len(old)} existing junction polygons...')
        commit = b.begin_commit()
        retry(b.remove_items, old)
        retry(b.push_commit, commit, 'Remove old LED junction polygons')
        shapes = retry(b.get_shapes)

    tracks = retry(b.get_tracks)
    vias = retry(b.get_vias)
    edge_stroke = max((s.attributes.stroke.width or 0 for s in shapes
                       if s.layer == bt.BoardLayer.Value('BL_Edge_Cuts')),
                      default=0)
    board_poly = board_area(b).buffer(-(EDGE_CLEARANCE + edge_stroke / 2),
                                      quad_segs=4)

    existing_poly_nets = {s.net.name for s in shapes
                          if isinstance(s, bt.BoardPolygon)
                          and s.layer in (F_CU, B_CU) and s.net.name}

    def on_layer(p, layer):
        return p['on_f'] if layer == F_CU else p['on_b']

    grow = CLEARANCE + CLEAR_EPS
    commit = b.begin_commit()
    made, new_items = [], []
    created = {F_CU: [], B_CU: []}
    try:
        for net_name in targets:
            if net_name not in junction_pads:
                print(f'SKIP {net_name}: no LED pads carry this net')
                continue
            if net_name in existing_poly_nets:
                print(f'SKIP {net_name}: polygon already exists')
                continue
            net = nets_by_name[net_name]
            group = junction_pads[net_name]
            refs = sorted({p['ref'] for p in group})

            pts = []
            for p in group:
                pts += list(p['box'].exterior.coords)

            end_union = {}
            for r in refs:
                end_union[r] = unary_union(
                    [p['box'] for p in group if p['ref'] == r])

            layer_polys, used_margin = None, None
            for margin in MARGINS:
                hull = octo_hull(pts, margin)
                # board edge: octilinear over-approximation of the forbidden
                # region (board_poly is already eroded by EDGE_CLEARANCE)
                edge_forbidden = hull.difference(board_poly)
                attempt = {}
                for layer in (F_CU, B_CU):
                    obstacles = []
                    if not edge_forbidden.is_empty:
                        obstacles.append(obstacle_oct(edge_forbidden, 0))
                    for p in pads:
                        if (p['eff'] == net_name or p['box'] is None
                                or not on_layer(p, layer)):
                            continue
                        if p['box'].distance(hull) < grow:
                            obstacles.append(obstacle_oct(p['box'], grow))
                    for t in tracks:
                        if (not isinstance(t, bt.Track) or t.layer != layer
                                or t.net.name == net_name):
                            continue
                        line = LineString([(t.start.x, t.start.y),
                                           (t.end.x, t.end.y)])
                        if line.distance(hull) < t.width / 2 + grow:
                            seg = line.intersection(
                                hull.buffer(t.width / 2 + grow, join_style=2))
                            if not seg.is_empty:
                                obstacles.append(obstacle_oct(
                                    seg.buffer(t.width / 2, quad_segs=8), grow))
                    for v in vias:
                        if v.net.name == net_name:
                            continue
                        g = Point(v.position.x, v.position.y)
                        if g.distance(hull) < v.diameter / 2 + grow:
                            obstacles.append(
                                obstacle_oct(g, v.diameter / 2 + grow))
                    for s in shapes:
                        if (isinstance(s, bt.BoardPolygon) and s.layer == layer
                                and s.net.name and s.net.name != net_name):
                            sgrow = grow + (s.attributes.stroke.width or 0) / 2
                            for pwh in s.polygons:
                                g = pwh_to_shapely(pwh)
                                if g.distance(hull) < sgrow:
                                    obstacles.append(obstacle_oct(g, sgrow))
                    for gnet, g in created[layer]:
                        if gnet != net_name and g.distance(hull) < grow:
                            obstacles.append(g.buffer(grow, join_style=2))

                    poly = (hull.difference(unary_union(obstacles))
                            if obstacles else hull)
                    pieces = (list(poly.geoms)
                              if isinstance(poly, MultiPolygon) else [poly])
                    keep = [pc for pc in pieces
                            if all(pc.intersects(end_union[r]) for r in refs)]
                    if not keep:
                        attempt = None
                        break
                    attempt[layer] = smooth(
                        max(keep, key=lambda pc: pc.area), end_union, refs)
                if attempt:
                    layer_polys, used_margin = attempt, margin
                    break

            if layer_polys is None:
                print(f'FAIL {net_name}: {refs[0]}<->{refs[1]} not bridgeable '
                      f'even with {to_mm(MARGINS[-1])} mm margin')
                continue

            for layer, poly in layer_polys.items():
                bp = bt.BoardPolygon()
                bp.layer = layer
                bp.net = net
                bp.polygons.append(shapely_to_pwh(poly))
                bp.attributes.fill.filled = True
                bp.attributes.stroke.width = 0
                new_items.append(bp)
                created[layer].append((net_name, poly))

            nv = [len(layer_polys[la].exterior.coords) - 1 for la in (F_CU, B_CU)]
            note = f' (margin {to_mm(used_margin):.1f})' if used_margin != MARGINS[0] else ''
            print(f'POLY {net_name}: {refs[0]}<->{refs[1]}  '
                  f'F.Cu {layer_polys[F_CU].area/1e12:.1f} mm2/{nv[0]}v, '
                  f'B.Cu {layer_polys[B_CU].area/1e12:.1f} mm2/{nv[1]}v{note}')
            made.append(net_name)

        if new_items:
            retry(b.create_items, new_items)
            retry(b.push_commit, commit,
                  f'LED junction copper polygons ({len(made)} nets)')
            print(f'\nCreated {len(new_items)} polygons for {len(made)} net(s).')
        else:
            b.drop_commit(commit)
            print('\nNothing created.')
    except Exception:
        b.drop_commit(commit)
        raise
    return 0


if __name__ == '__main__':
    sys.exit(main())
