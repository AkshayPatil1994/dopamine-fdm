#!/usr/bin/env python3
"""Recycled-precursor turbulence at the wave inlet (wave_turb = 1): a synthetic donor slice with spanwise-periodic fluctuations u', w' of amplitude 0.05
over the whole depth (the air included) is added to the wave + current inlet. The water receives it (none without wave_turb), the liquid ledger still
closes (1e-4). On the inlet face of the last snapshot, against the run without wave_turb at the same time: the injected u' = U_on - U_off has zero z-mean
(the donor mean is subtracted) and the rms and pattern of the donor slice in the water cells, the target w' = (W_ghost + W_first)/2 is the donor w' (peak
= amplitude) and both are exactly zero in the air above the free surface. max|W| of the monitor lies in (amp, 2 amp): the ghost mirror 2 w' - W_first, with
0 < W_first < amp because the synthetic slice is not divergence free and w' decays within a few cells."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='voftrb_')
n1, n1v, n2, nsn, amp = 32, 33, 5, 4, 0.05


kz = np.arange(n2) % (n2 - 1)   # the periodic last cell repeats the first
up = amp*np.cos(2*np.pi*(kz[:-1] + 0.5)/(n2 - 1))
wp = amp*np.sin(2*np.pi*kz[:-1]/(n2 - 1))


def donor(wd):
    u = np.tile(0.2 + amp*np.cos(2*np.pi*(kz + 0.5)/(n2 - 1)), (n1, 1))
    w = np.tile(amp*np.sin(2*np.pi*kz/(n2 - 1)), (n1, 1))
    v = np.zeros((n1v, n2))
    with open(os.path.join(wd, 'donor.bin'), 'wb') as f:
        for _ in range(nsn):
            for blk in (u, v, w):
                blk.T.copy().astype('>f8').tofile(f)   # big-endian (the build converts), j (first index) fastest
    np.linspace(0, 3, nsn).astype('>f8').tofile(os.path.join(wd, 'donor_times.bin'))
    open(os.path.join(wd, 'donor_meta.txt'), 'w').write('ncomp = 3\nn1 = %d\nn1_V = %d\nn2 = %d\ncomps = UVW\nnsnaps = %d\ndir = x\n' % (n1, n1v, n2, nsn))


def run(tag, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    text = re.sub(r'nsave = 100000', 'nsave = -1, tsave = 1.5', text, count=1)
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    donor(wd)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    d = np.loadtxt(os.path.join(wd, 'vof_diag.dat'))
    with open(os.path.join(wd, 'fields', 'flume.%d' % max(int(f.split('.')[1]) for f in os.listdir(os.path.join(wd, 'fields'))
                                                          if re.match(r'flume\.\d+$', f))), 'rb') as f:
        mesh = []
        for _ in range(6):
            n = int(np.fromfile(f, '>i4', 1)[0]); mesh.append(np.fromfile(f, '>f8', n))
        fld = []
        for _ in range(3):
            n = np.fromfile(f, '>i4', 3)
            fld.append(np.fromfile(f, '>f8', int(np.prod(n))).reshape(tuple(n), order='F'))
    return d[:, 17].max(), (d[-1, 3] - d[0, 3] - d[-1, 21] - d[-1, 22]) / d[0, 3], d[-1, 1], fld[0][0], fld[2][:2], mesh[4]


try:
    w_on, res, t_end, u_on, ww_on, ym = run('on', [])
    w_off, res_off, _, u_off, ww_off, _ = run('off', [(r'wave_turb = 1,', 'wave_turb = 0,')])
    nj, nk = len(ym), n2 - 1
    du = (u_on - u_off)[1:nj + 1, 1:nk + 1]            # cell j, z-cell k (first row / column: ghost)
    dw = (0.5*(ww_on[0] + ww_on[1]) - 0.5*(ww_off[0] + ww_off[1]))[1:nj + 1, :nk]   # w' target at the z-faces
    depth = 0.5
    wet, air = ym < depth - 0.02, ym > depth + 0.04    # cells wholly in the water / air at every phase of the wave
    e_mean = np.abs(du[wet].mean(axis=1)).max()
    e_u = np.abs(du[wet] - up).max()
    rms = np.sqrt((du[wet]**2).mean())
    e_w = np.abs(dw[wet] - wp).max()
    e_air = max(np.abs(du[air]).max(), np.abs(dw[air]).max())
    e_part = max(np.abs(du).max()/amp, np.abs(dw).max()/amp)
    ok = (t_end > 1.49 and abs(res) < 1e-4 and amp < w_on < 2*amp and w_off < 0.1*amp and e_mean < 1e-9 and e_u < 1e-6
          and abs(rms/np.sqrt((up**2).mean()) - 1) < 1e-6 and e_w < 1e-6 and e_air < 1e-12 and e_part <= 1 + 1e-6)
    print('max|W| with %.4f without %.4f (amplitude %.2f)  ledger residual %.2e (%.2e without)' % (w_on, w_off, amp, res, res_off))
    print("inlet face: z-mean of u' %.1e  |u' - donor| %.1e  rms(u')/rms(donor) %.6f  |w' - donor| %.1e  in the air %.1e  %s"
          % (e_mean, e_u, rms/np.sqrt((up**2).mean()), e_w, e_air, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
