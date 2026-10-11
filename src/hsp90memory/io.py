"""Loading and idealising plasmon-ruler traces (Zenodo 12697153 format: one column of raw counts, 10 Hz)."""
import numpy as np

DT = 0.1  # s, frame interval


def load_trace(path):
    """Raw counts -> relative intensity I/<I>."""
    raw = np.loadtxt(path)
    return raw / raw.mean()


def idealise(I):
    """Two-level threshold at the midpoint between the open and closed peaks.
    At SNR ~ 11 the false-crossing probability is ~1e-8 per frame, so this is error-free for these data.
    Returns an int array, 0 = open (low intensity), 1 = closed (high)."""
    thr = (np.percentile(I, 20) + np.percentile(I, 95)) / 2
    return (I > thr).astype(np.int8)


def dwells(s):
    """Dwell lengths (frames) and classes of a 0/1 series, including the censored first and last dwells."""
    j = np.flatnonzero(np.diff(s)) + 1
    b = np.r_[0, j, s.size]
    return np.diff(b), s[b[:-1]]


def load_dwells(path):
    s = idealise(load_trace(path))
    dw, dc = dwells(s)
    return s, dw, dc
