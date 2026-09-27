"""Helper run with the system python3 (which has pcbnew): dump exact pad
geometry from a board file as JSON.  Works around a KiCad IPC API bug where
get_pad_shapes_as_polygons returns offset shapes for rotated footprints.

Usage: python3 pad_dump_helper.py BOARD_FILE > pads.json
"""
import json
import sys

import pcbnew

board = pcbnew.LoadBoard(sys.argv[1])
out = []
for fp in board.GetFootprints():
    ref = fp.GetReference()
    fpid = str(fp.GetFPID().GetLibItemName())
    for p in fp.Pads():
        try:
            poly = p.GetEffectivePolygon()
        except TypeError:
            poly = p.GetEffectivePolygon(pcbnew.ERROR_INSIDE)
        pts = []
        if poly.OutlineCount() > 0:
            o = poly.Outline(0)
            pts = [[o.CPoint(i).x, o.CPoint(i).y] for i in range(o.PointCount())]
        out.append({
            'ref': ref,
            'fpid': fpid,
            'num': str(p.GetNumber()),
            'net': p.GetNetname(),
            'on_f': p.IsOnLayer(pcbnew.F_Cu),
            'on_b': p.IsOnLayer(pcbnew.B_Cu),
            'poly': pts,
        })
json.dump(out, sys.stdout)
