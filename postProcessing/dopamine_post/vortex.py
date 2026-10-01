"""
vortex.py -- Q-criterion vortex identification and iso-surface movies.

  ANALYSIS   q_criterion(snap)        Q = 0.5 (|Omega|^2 - |S|^2) of a FieldSnapshot's
                                      cell-centred velocity (non-uniform grid ok).
             q_volume_fraction(...)   fraction of the domain with Q > threshold.
  MOVIE      animate_q_isosurface     Q iso-surface coloured by |u| next to z- and
                                      x-normal |u| slices (animateQ.py).

The movie needs VTK (`pip install vtk`) and an X display (use `xvfb-run` on a
headless node); ffmpeg must be on PATH.  The analysis half needs only numpy.

Library use
-----------
    from dopamine_post.fields import FieldSeries
    from dopamine_post.vortex import q_criterion, animate_q_isosurface

    series = FieldSeries(".")
    Q = q_criterion(series.latest())
    animate_q_isosurface(series, "q.mp4", qval=5.0, start=1000, stride=2)
"""
import subprocess

import numpy as np

from . import _core

# camera direction (focal point -> camera, in units of Lx) and view-up
VIEWS = {
    "iso":   ((-0.75, 0.55, 0.85), (0, 1, 0)),
    "side":  ((0, 0, 1), (0, 1, 0)),        # x-y plane, looking along -z
    "top":   ((0, 1, 0), (0, 0, -1)),       # x-z plane, looking down -y
    "front": ((-1, 0, 0), (0, 1, 0)),       # end-on, looking along +x
}

_LAYOUTS = {
    "both":  [("q", (0.0, 0, 0.58, 1)), ("zslice", (0.58, 0.80, 1, 1)), ("xslice", (0.58, 0, 1, 0.80))],
    "q":     [("q", (0, 0, 1, 1))],
    "slice": [("zslice", (0, 0.80, 1, 1)), ("xslice", (0, 0, 1, 0.80))],
}


# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS
# ════════════════════════════════════════════════════════════════════════════
def q_criterion(snap):
    """Q = 0.5 (|Omega|^2 - |S|^2) on the cell centres of a FieldSnapshot, shape (nxm, nym, nzm).

    Gradients are second-order central differences on the (possibly stretched) cell-centre
    grid, one-sided at the domain edges.  Q > 0 marks rotation-dominated regions.
    """
    vel = (snap["U"], snap["V"], snap["W"])
    coords = (snap.xm, snap.ym, snap.zm)
    g = [[np.gradient(vel[i], coords[j], axis=j) for j in range(3)] for i in range(3)]
    Q = np.zeros_like(vel[0])
    for i in range(3):
        for j in range(3):
            S = 0.5 * (g[i][j] + g[j][i])
            O = 0.5 * (g[i][j] - g[j][i])
            Q += 0.5 * (O * O - S * S)
    return Q


def q_volume_fraction(snap, qval):
    """Volume fraction of the domain where Q > qval (cell volumes from the face grid)."""
    Q = q_criterion(snap)
    dV = (np.diff(snap.x)[:, None, None] * np.diff(snap.y)[None, :, None]
          * np.diff(snap.z)[None, None, :])
    return float((dV * (Q > qval)).sum() / dV.sum())


# ════════════════════════════════════════════════════════════════════════════
# MOVIE
# ════════════════════════════════════════════════════════════════════════════
def _import_vtk():
    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError as e:
        raise ImportError("animate_q_isosurface needs VTK: pip install vtk") from e
    return vtk, numpy_support


def _make_grid(vtk, ns, xm, ym, zm, arrays):
    grid = vtk.vtkRectilinearGrid()
    grid.SetDimensions(len(xm), len(ym), len(zm))
    for setter, c in ((grid.SetXCoordinates, xm), (grid.SetYCoordinates, ym),
                      (grid.SetZCoordinates, zm)):
        setter(ns.numpy_to_vtk(np.ascontiguousarray(c, dtype=np.float64)))
    for name, arr in arrays.items():
        # VTK point order is x fastest -> Fortran-order flatten of (nx, ny, nz)
        v = ns.numpy_to_vtk(arr.ravel(order="F").astype(np.float32), deep=True)
        v.SetName(name)
        grid.GetPointData().AddArray(v)
    grid.GetPointData().SetActiveScalars("Q")
    return grid


def _make_lut(vtk, vmin, vmax, cmap):
    from matplotlib import colormaps
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetRange(vmin, vmax)
    cm = colormaps[cmap]
    for i in range(256):
        lut.SetTableValue(i, *cm(i / 255.0))
    lut.Build()
    return lut


def _outline(vtk, grid):
    f = vtk.vtkOutlineFilter()
    f.SetInputData(grid)
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(f.GetOutputPort())
    a = vtk.vtkActor()
    a.SetMapper(m)
    a.GetProperty().SetColor(0, 0, 0)
    return a


def _colour_mapper(vtk, port, lut, name):
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(port)
    m.SetLookupTable(lut)
    m.SetScalarRange(*lut.GetRange())
    m.SelectColorArray(name)
    m.SetScalarModeToUsePointFieldData()
    return m


def _add_isosurface(vtk, ren, grid, lut, qval):
    c = vtk.vtkContourFilter()
    c.SetInputData(grid)
    c.SetInputArrayToProcess(0, 0, 0, 0, "Q")
    c.SetValue(0, qval)
    c.ComputeNormalsOn()
    c.ComputeScalarsOff()
    probe = vtk.vtkProbeFilter()            # carry |u| onto the surface
    probe.SetInputConnection(c.GetOutputPort())
    probe.SetSourceData(grid)
    a = vtk.vtkActor()
    a.SetMapper(_colour_mapper(vtk, probe.GetOutputPort(), lut, "Umag"))
    a.GetProperty().SetSpecular(0.3)
    ren.AddActor(a)


def _add_slice(vtk, ren, grid, lut, origin, normal):
    pl = vtk.vtkPlane()
    pl.SetOrigin(*origin)
    pl.SetNormal(*normal)
    cut = vtk.vtkCutter()
    cut.SetCutFunction(pl)
    cut.SetInputData(grid)
    a = vtk.vtkActor()
    a.SetMapper(_colour_mapper(vtk, cut.GetOutputPort(), lut, "Umag"))
    a.GetProperty().LightingOff()
    ren.AddActor(a)


def _add_axes(vtk, rw, viewport, main_ren):
    ren = vtk.vtkRenderer()
    ren.SetViewport(*viewport)
    ren.SetLayer(1)
    ren.InteractiveOff()
    ren.SetActiveCamera(main_ren.GetActiveCamera())
    ax = vtk.vtkAxesActor()
    ax.SetTotalLength(1, 1, 1)
    ren.AddActor(ax)
    rw.AddRenderer(ren)


def _title(vtk, text):
    t = vtk.vtkTextActor()
    t.SetInput(text)
    t.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    t.GetPositionCoordinate().SetValue(0.01, 0.98)
    t.GetTextProperty().SetVerticalJustificationToTop()
    t.GetTextProperty().SetColor(0.1, 0.1, 0.1)
    t.GetTextProperty().SetFontSize(14)
    return t


def _setup_camera(ren, centre, L, view, cam_pos, azimuth, elevation, roll, zoom, parallel):
    direction, up = VIEWS[view]
    d = np.array(cam_pos if cam_pos else direction, float)
    d /= np.linalg.norm(d)
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(*centre)
    cam.SetViewUp(*up)
    cam.SetPosition(*(np.array(centre) + d * L))
    if parallel:
        cam.ParallelProjectionOn()
    ren.ResetCamera()
    cam.Azimuth(azimuth)
    cam.Elevation(elevation)
    cam.OrthogonalizeViewUp()
    cam.Roll(roll)
    cam.Zoom(zoom)


def _setup_slice_camera(ren, centre, normal, width, height, aspect, zoom):
    """Head-on orthographic camera fitting a width x height plane, leaving room for the title."""
    cam = ren.GetActiveCamera()
    vis_h = max(height * 1.14, width * 1.03 / aspect) / zoom
    shift = (vis_h - height) / 2 - 0.01 * vis_h
    centre = np.array(centre, float) + np.array([0, shift, 0])
    cam.ParallelProjectionOn()
    cam.SetFocalPoint(*centre)
    cam.SetPosition(*(centre + 10 * max(width, height) * np.array(normal, float)))
    cam.SetViewUp(0, 1, 0)
    ren.ResetCameraClippingRange()
    cam.SetParallelScale(vis_h / 2)


def animate_q_isosurface(series, out="q.mp4", steps=None, start=None, end=None, stride=1,
                         max_frames=None, qval=5.0, vmax=None, cmap="RdYlBu_r", fps=10,
                         width=None, height=600, plot="both", view="iso", cam_pos=None,
                         azimuth=0.0, elevation=0.0, roll=0.0, zoom=1.3, zoom_slice=1.0,
                         slice_x=0.5, parallel=False):
    """Animate the Q = `qval` iso-surface (coloured by |u|) of a FieldSeries (animateQ.py).

    plot='both' adds a z-normal (mid-span) and an x-normal (at `slice_x` of Lx) |u| slice
    beside the iso-surface; 'q' / 'slice' draws only one of the two.  Snapshots are chosen
    either from an explicit `steps` list or by `start`/`end` step bounds and `stride`.
    `vmax` is the top of the |u| colour range (default 1.7 * Ub_target).  `view` is one of
    VIEWS ('iso', 'side', 'top', 'front'); `cam_pos` (dx, dy, dz) overrides it.  Output is
    .mp4 or .gif, chosen by the extension of `out`; rendering is off-screen VTK, so it
    still needs an X display.
    """
    vtk, ns = _import_vtk()
    out = str(out)
    snaps = series.snapshots
    if steps is not None:
        keep = set(steps)
        snaps = [s for s in snaps if s[0] in keep]
    else:
        snaps = [s for s in snaps if (start is None or s[0] >= start) and (end is None or s[0] <= end)]
        snaps = snaps[::max(1, stride)]
    if max_frames:
        snaps = snaps[:max_frames]
    if not snaps:
        raise FileNotFoundError(f"no snapshots selected in {series.fields_dir}")
    print(f"{len(snaps)} frames -> {out}")

    ub = float(_core.read_input(series.case_dir, ("Ub_target",)).get("Ub_target", 1.0))
    lut = _make_lut(vtk, 0.0, vmax or 1.7 * ub, cmap)

    panels = _LAYOUTS[plot]
    width = width or (1360 if plot == "both" else 680)
    w, h = width - width % 2, height - height % 2
    rw = vtk.vtkRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetSize(width, height)
    rw.SetMultiSamples(4)
    rw.SetNumberOfLayers(2)
    titles = {"q": f"Q = {qval:g} iso-surface, coloured by velocity magnitude",
              "zslice": "z-normal slice (mid-span)", "xslice": "x-normal slice"}
    rens = []
    for kind, vp in panels:
        r = vtk.vtkRenderer()
        r.SetViewport(*vp)
        r.SetBackground(0.99, 0.99, 0.98)
        r.AddViewProp(_title(vtk, titles[kind]))
        rw.AddRenderer(r)
        rens.append(r)

    vf = "vflip" if out.endswith(".mp4") else "vflip,split[a][b];[a]palettegen[p];[b][p]paletteuse"
    codec = ["-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18"] if out.endswith(".mp4") else ["-loop", "0"]
    ffmpeg = subprocess.Popen(
        ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
         "-s", f"{w}x{h}", "-r", str(fps), "-i", "-", "-vf", vf] + codec + [out],
        stdin=subprocess.PIPE)
    grab = vtk.vtkWindowToImageFilter()
    grab.SetInput(rw)
    grab.ReadFrontBufferOff()

    from .fields import FieldSnapshot
    grid = None
    for k, (step, path) in enumerate(snaps):
        snap = FieldSnapshot.read(path)
        xm, ym, zm = snap.xm, snap.ym, snap.zm
        mag = np.sqrt(snap["U"] ** 2 + snap["V"] ** 2 + snap["W"] ** 2)
        g = _make_grid(vtk, ns, xm, ym, zm, {"Q": q_criterion(snap), "Umag": mag})
        if grid is None:
            grid = g
            x0 = xm[0] + slice_x * (xm[-1] - xm[0])
            z0 = 0.5 * (zm[0] + zm[-1])
            centre = (xm.mean(), ym.mean(), zm.mean())
            L = xm[-1] - xm[0]
            for (kind, vp), r in zip(panels, rens):
                r.AddActor(_outline(vtk, grid))
                pw, ph = (vp[2] - vp[0]) * width, (vp[3] - vp[1]) * height
                if kind == "q":
                    _add_isosurface(vtk, r, grid, lut, qval)
                    _setup_camera(r, centre, L, view, cam_pos, azimuth, elevation, roll, zoom, parallel)
                    _add_axes(vtk, rw, (vp[0], vp[1], vp[0] + 0.08, vp[1] + 0.16), r)
                elif kind == "zslice":
                    _add_slice(vtk, r, grid, lut, (0, 0, z0), (0, 0, 1))
                    _setup_slice_camera(r, (centre[0], centre[1], z0), (0, 0, 1),
                                        L, ym[-1] - ym[0], pw / ph, zoom_slice)
                else:
                    _add_slice(vtk, r, grid, lut, (x0, 0, 0), (1, 0, 0))
                    _setup_slice_camera(r, (x0, centre[1], centre[2]), (-1, 0, 0),
                                        zm[-1] - zm[0], ym[-1] - ym[0], pw / ph, zoom_slice)
        else:
            grid.DeepCopy(g)
        grid.Modified()

        rw.Render()
        grab.Modified()
        grab.Update()
        img = grab.GetOutput()
        dims = img.GetDimensions()
        arr = ns.vtk_to_numpy(img.GetPointData().GetScalars()).reshape(dims[1], dims[0], -1)[:h, :w, :3]
        ffmpeg.stdin.write(np.ascontiguousarray(arr).tobytes())
        print(f"  [{k + 1}/{len(snaps)}] step {step}", flush=True)

    ffmpeg.stdin.close()
    ffmpeg.wait()
    print(f"Wrote {out}")
