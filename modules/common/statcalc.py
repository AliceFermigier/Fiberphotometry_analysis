# -*- coding: utf-8 -*-

import numpy as np

def epochs_from_indices(sig, idx, pre, post, sr, bl=(2.0, 0.0), std0=None):
    """Return (n_events, T) baseline-subtracted epochs. std0=None -> raw dFF; std0=float -> z-score."""
    T = int(round((pre + post) * sr)) + 1
    out = np.full((len(idx), T), np.nan)
    for i, e in enumerate(idx):
        s, t = int(e - pre * sr), int(e + post * sr)
        b0, b1 = int(e - bl[0] * sr), int(e - bl[1] * sr)
        if s < 0 or t >= len(sig) or b0 < 0:
            continue
        F0 = np.nanmean(sig[b0:b1 + 1])
        ep = sig[s:t + 1] - F0
        out[i] = ep / std0 if std0 else ep
    return out


def shuffle_null(sig, event_idx, allowed_mask, pre, post, sr, bl=(2.0, 0.0),
                 std0=None, n_perm=2000, seed=0):
    """
    Mean PETH under random event times.
    allowed_mask: bool array of frames where a random 'event' may be placed
                  (e.g. valid & no scored behaviour).
    """
    rng = np.random.default_rng(seed)
    pad = int(max(pre, bl[0]) * sr), int(post * sr)
    cand = np.where(allowed_mask)[0]
    cand = cand[(cand >= pad[0]) & (cand < len(sig) - pad[1])]
    n_ev = len(event_idx)
    null = np.empty((n_perm, int(round((pre + post) * sr)) + 1))
    for k in range(n_perm):
        rnd = rng.choice(cand, n_ev, replace=False)
        null[k] = np.nanmean(epochs_from_indices(sig, rnd, pre, post, sr, bl, std0), axis=0)
    return null


def significance_from_null(obs_mean, null, alpha=0.05, min_consec=3, sr=10, post_mask=None):
    """
    Pointwise one-sided test (rise above chance) + family-wise max-statistic threshold.
    Returns dict of boolean arrays: 'pointwise', 'consecutive', 'familywise'.
    """
    if post_mask is None:
        post_mask = np.ones(obs_mean.size, dtype=bool)
    thr_pt = np.nanpercentile(null, 100 * (1 - alpha), axis=0)
    thr_fw = np.nanpercentile(np.nanmax(null[:, post_mask], axis=1), 100 * (1 - alpha))

    pointwise = (obs_mean > thr_pt) & post_mask
    # require min_consec successive significant samples
    consec = np.zeros_like(pointwise)
    run = 0
    for i, v in enumerate(pointwise):
        run = run + 1 if v else 0
        if run >= min_consec:
            consec[i - min_consec + 1:i + 1] = True
    familywise = (obs_mean > thr_fw) & post_mask
    return {'pointwise': pointwise, 'consecutive': consec, 'familywise': familywise,
            'thr_pt': thr_pt, 'thr_fw': thr_fw}
