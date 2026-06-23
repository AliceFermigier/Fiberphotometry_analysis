import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import modules.common.preprocess as pp

def txt_to_df(capacitance_txt_path):
    # Load all columns
    df = pd.read_csv(
        capacitance_txt_path,
        header=None,
        names=["time(ms)", "Capacitance", "recording"],
        low_memory=False
    )

    # Find indices for recording start and stop
    try:
        start_idx = df.index[df["recording"] == "RECORDING_START"][0]
        stop_idx = df.index[df["recording"] == "RECORDING_STOP"][0]
    except IndexError:
        raise ValueError("RECORDING_START or RECORDING_STOP not found in file")

    # Slice the dataframe to include only the recording period
    df = df.loc[start_idx:stop_idx]

    # Reset time to start at 0 and convert to seconds
    df["Time(s)"] = (df["time(ms)"] - df["time(ms)"].iloc[0]) / 1000

    # Drop the original millisecond column
    df = df[["Time(s)", "Capacitance"]]

    # Convert capacitance to integer
    df["Capacitance"] = df["Capacitance"].astype(int)

    return df

def extract_lick_bouts(capacitance_df, threshold):
    licks_df = (
        capacitance_df[capacitance_df["Capacitance"] > threshold][["Time(s)"]]
        .reset_index(drop=True)
    )
    return licks_df

def plot_licks_and_threshold(capacitance_df, threshold=500, downsample=100):

    capacitance_df_downsampled = pp.downsample(capacitance_df, target_frequency=downsample)

    fig = px.line(
        capacitance_df_downsampled,
        x="Time(s)",
        y="Capacitance",
        title="Capacitance over time",
        labels={"Time(s)": "Time (s)", "Capacitance": "Capacitance"},
    )

    # Threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="red",
        annotation_text=f"Threshold ({threshold})",
        annotation_position="top right",
    )

    fig.update_layout(
        xaxis_title="Time (s)",
        yaxis_title="Capacitance",
        hovermode="x unified",   # shows all values at a given x on hover
    )

    fig.show()

def compute_distance(df, port, scale, nose_x='nose_x', nose_y='nose_y'):
    """
    Compute distance from nose to port in centimeters.
    
    Parameters:
    -----------
    df : DataFrame
        Contains nose position columns (in pixels)
    port : dict
        Port coordinates with 'x' and 'y' keys (in pixels)
    scale : float
        Scale factor from scale_and_coords["Scale_cm_per_px"]
    nose_x, nose_y : str
        Column names for nose coordinates (in pixels)
    
    Returns:
    --------
    Distance in centimeters
    """
    px = port["x"]
    py = port["y"]
    pixel_dist = np.sqrt((df[nose_x] - px)**2 + (df[nose_y] - py)**2)
    return pixel_dist * scale  # Convert pixels to cm

def filter_licking(df, ports, scale_and_coords, lick_col='Licks', lick_radius_cm=0.5):
    """
    Filter licks based on distance from lick port.
    
    Parameters:
    -----------
    lick_radius_cm : float
        Radius in centimeters within which a lick is valid 
    """
    video_scale = scale_and_coords["Scale_cm_per_px"]
    dist = compute_distance(df, ports["lick_port"], video_scale)
    true_lick = (df[lick_col] == 1) & (dist < lick_radius_cm)
    df["Licks_filtered"] = true_lick.astype(int)
    return df

def detect_airpuff_entry(df, ports, scale_and_coords, radius_cm=3.0):
    """
    Detect when nose enters airpuff ports.
    
    Parameters:
    -----------
    radius_cm : float
        Radius in centimeters defining port entry zone
    """
    video_scale = scale_and_coords["Scale_cm_per_px"]
    distL = compute_distance(df, ports["airpuff_left"], video_scale)
    distR = compute_distance(df, ports["airpuff_right"], video_scale)
    df["Nose_in_any_airport"] = ((distL < radius_cm) | (distR < radius_cm)).astype(int)
    return df

def detect_approach(df, ports, scale_and_coords,
                    approach_radius_cm=10.0,
                    speed_col='Speed'):
    """
    Detect when the mouse is approaching a port: within radius AND moving
    toward it (distance decreasing over time).

    A frame is flagged as 'approaching' if:
      1. Nose is within approach_radius_cm of the port.
      2. The derivative of distance is negative (closing in).

    Parameters
    ----------
    df : pd.DataFrame
        Must contain nose coordinates, Time(s), and optionally Speed.
    ports : dict
        Port coordinates dict with keys 'airpuff_left', 'airpuff_right', 'lick_port'.
    scale_and_coords : dict
        Must contain 'Scale_cm_per_px'.
    approach_radius_cm : float
        Outer radius defining the approach zone (cm).
    speed_col : str
        Speed column name — used to exclude immobile frames (speed = 0 is
        not an approach, it's a stop).

    Returns
    -------
    df with new columns:
        'Approach_airpuff_left'   — 1 when actively approaching left port
        'Approach_airpuff_right'  — 1 when actively approaching right port
        'Approach_any_airpuff'    — 1 when approaching either airpuff port
        'Approach_lickport'       — 1 when actively approaching lick port
    """
    scale = scale_and_coords["Scale_cm_per_px"]

    port_map = {
        'airpuff_left'  : 'Approach_airpuff_left',
        'airpuff_right' : 'Approach_airpuff_right',
        'lick_port'     : 'Approach_lickport',
    }

    for port_key, col_name in port_map.items():
        dist = compute_distance(df, ports[port_key], scale)

        # Derivative of distance: negative = moving closer
        dist_diff = dist.diff().fillna(0)

        within_radius = dist < approach_radius_cm
        moving_closer = dist_diff < 0

        # Exclude frames where mouse is stopped (speed = 0)
        if speed_col in df.columns:
            moving = df[speed_col] > 0
        else:
            moving = pd.Series(True, index=df.index)

        df[col_name] = (within_radius & moving_closer & moving).astype(int)

    df['Approach_any_airpuff'] = (
        (df['Approach_airpuff_left'] == 1) | (df['Approach_airpuff_right'] == 1)
    ).astype(int)

    return df


def detect_head_orientation(df, ports, scale_and_coords,
                             angle_threshold_deg=45.0,
                             nose_x='nose_x', nose_y='nose_y',
                             ear_left_x='ear_left_x', ear_left_y='ear_left_y',
                             ear_right_x='ear_right_x', ear_right_y='ear_right_y'):
    """
    Detect when the mouse's head is directed toward a port.

    Head direction is the vector from the midpoint between the two ears
    (proxy for the back of the head) to the nose. The angle between this
    vector and the nose-to-port vector is computed. If the angle is below
    angle_threshold_deg the head is considered oriented toward the port.

    Parameters
    ----------
    angle_threshold_deg : float
        Maximum angle (degrees) between head direction and port direction
        for the head to be considered 'oriented' toward the port.

    Returns
    -------
    df with new columns:
        'Head_toward_airpuff_left'
        'Head_toward_airpuff_right'
        'Head_toward_any_airpuff'
        'Head_toward_lickport'
        'Angle_to_airpuff_left'    — raw angle in degrees (useful for QC)
        'Angle_to_airpuff_right'
        'Angle_to_lickport'
    """
    # Ear midpoint = proxy for back of head
    ear_mid_x = (df[ear_left_x] + df[ear_right_x]) / 2
    ear_mid_y = (df[ear_left_y] + df[ear_right_y]) / 2

    # Head direction vector: ear midpoint → nose
    head_vec_x = df[nose_x] - ear_mid_x
    head_vec_y = df[nose_y] - ear_mid_y

    def _angle_to_port(port):
        # Vector from nose to port
        port_vec_x = port['x'] - df[nose_x]
        port_vec_y = port['y'] - df[nose_y]

        # Dot product and magnitudes
        dot = head_vec_x * port_vec_x + head_vec_y * port_vec_y
        mag_head = np.sqrt(head_vec_x**2 + head_vec_y**2)
        mag_port = np.sqrt(port_vec_x**2 + port_vec_y**2)

        # Clamp to [-1, 1] to guard against floating point errors in arccos
        cos_angle = np.clip(dot / (mag_head * mag_port + 1e-8), -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))

    port_map = {
        'airpuff_left'  : ('Angle_to_airpuff_left',  'Head_toward_airpuff_left'),
        'airpuff_right' : ('Angle_to_airpuff_right', 'Head_toward_airpuff_right'),
        'lick_port'     : ('Angle_to_lickport',      'Head_toward_lickport'),
    }

    for port_key, (angle_col, oriented_col) in port_map.items():
        angles = _angle_to_port(ports[port_key])
        df[angle_col]    = angles
        df[oriented_col] = (angles < angle_threshold_deg).astype(int)

    df['Head_toward_any_airpuff'] = (
        (df['Head_toward_airpuff_left'] == 1) | (df['Head_toward_airpuff_right'] == 1)
    ).astype(int)

    return df

def detect_stretch(df, ports, scale_and_coords,
                   nose_center_threshold_cm=3.0,
                   orientation_angle_threshold_deg=60.0,
                   nose_x='nose_x', nose_y='nose_y',
                   center_x='center_x', center_y='center_y',
                   ear_left_x='ear_left_x', ear_left_y='ear_left_y',
                   ear_right_x='ear_right_x', ear_right_y='ear_right_y'):
    """
    Detect stretching posture toward airpuff ports.

    Criteria (both must be true):
      1. EXTENSION — nose-to-center distance > nose_center_threshold_cm,
         indicating the mouse is reaching its nose forward.
      2. ORIENTATION — head vector is pointed toward the port within
         orientation_angle_threshold_deg.

    Parameters
    ----------
    nose_center_threshold_cm : float
        Minimum nose-to-center distance (cm) to flag as extended/stretching.
        Plot df['Nose_center_dist_cm'].hist() on a few sessions to tune.
    """
    scale = scale_and_coords["Scale_cm_per_px"]

    # ── Nose-to-center distance ───────────────────────────────────────────────
    nose_center_dist_px = np.sqrt(
        (df[nose_x] - df[center_x])**2 +
        (df[nose_y] - df[center_y])**2
    )
    nose_center_dist_cm = nose_center_dist_px * scale
    df['Nose_center_dist_cm'] = nose_center_dist_cm          # keep for QC/tuning
    extended = nose_center_dist_cm > nose_center_threshold_cm

    # ── Head orientation ──────────────────────────────────────────────────────
    ear_mid_x  = (df[ear_left_x] + df[ear_right_x]) / 2
    ear_mid_y  = (df[ear_left_y] + df[ear_right_y]) / 2
    head_vec_x = df[nose_x] - ear_mid_x
    head_vec_y = df[nose_y] - ear_mid_y

    def _oriented_toward(port):
        port_vec_x = port['x'] - df[nose_x]
        port_vec_y = port['y'] - df[nose_y]
        dot        = head_vec_x * port_vec_x + head_vec_y * port_vec_y
        mag        = (np.sqrt(head_vec_x**2 + head_vec_y**2) *
                      np.sqrt(port_vec_x**2 + port_vec_y**2) + 1e-8)
        return np.degrees(np.arccos(np.clip(dot / mag, -1.0, 1.0))) < orientation_angle_threshold_deg

    for port_key, col in [('airpuff_left',  'Stretch_toward_airpuff_left'),
                           ('airpuff_right', 'Stretch_toward_airpuff_right')]:
        df[col] = (extended & _oriented_toward(ports[port_key])).astype(int)

    df['Stretch_toward_any_airpuff'] = (
        (df['Stretch_toward_airpuff_left'] == 1) | (df['Stretch_toward_airpuff_right'] == 1)
    ).astype(int)

    return df

def detect_rearing_and_grooming(df, scale_and_coords,
                                 body_compression_threshold_cm,
                                 nose_jitter_threshold_cm=0.3,
                                 ear_spread_threshold_cm=None,
                                 jitter_window_frames=10,
                                 speed_threshold_cms=2.0,
                                 nose_x='nose_x', nose_y='nose_y',
                                 tail_base_x='tail_base_x', tail_base_y='tail_base_y',
                                 center_x='center_x', center_y='center_y',
                                 ear_left_x='ear_left_x', ear_left_y='ear_left_y',
                                 ear_right_x='ear_right_x', ear_right_y='ear_right_y',
                                 speed_col='Speed'):
    """
    Detect rearing and grooming from a TOP-DOWN camera.

    Both behaviours share body compression (low nose-tailbase distance) and
    low locomotion speed. They are separated by nose micro-movement and ear
    spread:

      - REARING   : compressed + slow + LOW nose jitter + stable ear spread
      - GROOMING  : compressed + slow + HIGH nose jitter + narrow ear spread

    !! IMPORTANT — thresholds are dataset-dependent !!
    Run once with body_compression_threshold_cm set and
    nose_jitter_threshold_cm / ear_spread_threshold_cm as None to print
    data-driven suggestions, then validate against hand-scored video.

    Parameters
    ----------
    body_compression_threshold_cm : float
        Nose-to-tailbase distance (cm) below which body is considered compressed.
        Shared gate for both rearing and grooming detection.
    nose_jitter_threshold_cm : float
        Rolling std of nose position (cm) separating grooming (high) from
        rearing (low). Tune from df['Nose_jitter_cm'].hist().
    ear_spread_threshold_cm : float or None
        Inter-ear distance (cm) below which ears are considered narrow
        (consistent with facial grooming). If None, uses the 25th percentile
        of the session as a suggestion.
    jitter_window_frames : int
        Rolling window (frames) for computing nose jitter.
        Should span roughly 0.5 s (e.g. 10 frames at 20 fps).
    speed_threshold_cms : float
        Speed below which the mouse is considered stationary.

    Returns
    -------
    df with new columns:
        'Nose_tailbase_dist_cm'  — primary compression metric
        'Nose_jitter_cm'         — local nose movement (rolling std)
        'Ear_spread_cm'          — inter-ear distance
        'Body_compressed'        — gate: 1 when body is compressed + slow
        'Rearing'                — 1 when rearing
        'Grooming'               — 1 when grooming
    """
    scale = scale_and_coords["Scale_cm_per_px"]

    # ── Nose-to-tailbase distance (shared compression gate) ───────────────────
    nose_tail_px = np.sqrt(
        (df[nose_x] - df[tail_base_x])**2 +
        (df[nose_y] - df[tail_base_y])**2
    )
    df['Nose_tailbase_dist_cm'] = nose_tail_px * scale
    compressed = df['Nose_tailbase_dist_cm'] < body_compression_threshold_cm

    # ── Speed gate ────────────────────────────────────────────────────────────
    if speed_col in df.columns:
        slow = df[speed_col] < speed_threshold_cms
    else:
        print("  [detect] Speed column not found — omitting speed gate.")
        slow = pd.Series(True, index=df.index)

    body_compressed = compressed & slow
    df['Body_compressed'] = body_compressed.astype(int)

    # ── Nose jitter: rolling std of nose position ─────────────────────────────
    nose_pos_cm = pd.Series(
        np.sqrt(df[nose_x]**2 + df[nose_y]**2) * scale,
        index=df.index
    )
    df['Nose_jitter_cm'] = (
        nose_pos_cm
        .rolling(jitter_window_frames, center=True, min_periods=1)
        .std()
        .fillna(0)
    )

    if nose_jitter_threshold_cm is None:
        suggested = np.nanpercentile(
            df.loc[body_compressed, 'Nose_jitter_cm'], 50)
        print(f"  [detect] nose_jitter_threshold_cm not set. "
              f"Median jitter during compressed frames = {suggested:.3f} cm. "
              f"Plot df['Nose_jitter_cm'].hist() to tune.")
        nose_jitter_threshold_cm = suggested

    high_jitter = df['Nose_jitter_cm'] >= nose_jitter_threshold_cm
    low_jitter  = ~high_jitter

    # ── Ear spread: inter-ear distance ────────────────────────────────────────
    ear_spread_px = np.sqrt(
        (df[ear_left_x] - df[ear_right_x])**2 +
        (df[ear_left_y] - df[ear_right_y])**2
    )
    df['Ear_spread_cm'] = ear_spread_px * scale

    if ear_spread_threshold_cm is None:
        suggested = np.nanpercentile(
            df.loc[body_compressed, 'Ear_spread_cm'], 25)
        print(f"  [detect] ear_spread_threshold_cm not set. "
              f"25th percentile of ear spread during compressed frames = "
              f"{suggested:.2f} cm. Plot df['Ear_spread_cm'].hist() to tune.")
        ear_spread_threshold_cm = suggested

    narrow_ears = df['Ear_spread_cm'] < ear_spread_threshold_cm

    # ── Classification ────────────────────────────────────────────────────────
    # Rearing  : compressed + low jitter (head still)
    # Grooming : compressed + high jitter + narrow ears
    # Ambiguous compressed frames (high jitter but normal ears) are left as 0
    df['Rearing']  = (body_compressed & low_jitter).astype(int)
    df['Grooming'] = (body_compressed & high_jitter & narrow_ears).astype(int)

    # Rearing takes priority if both fire (shouldn't happen often)
    df.loc[df['Rearing'] == 1, 'Grooming'] = 0

    for label in ['Rearing', 'Grooming']:
        n = df[label].sum()
        print(f"  [detect] {label}: {n} frames ({n / len(df) * 100:.1f}%)")

    return df