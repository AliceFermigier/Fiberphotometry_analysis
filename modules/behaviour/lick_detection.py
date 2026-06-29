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
                    arrival_threshold_cm=1.0,
                    min_approach_duration_s=1.0,
                    max_distance_increase_cm=2.0,
                    smoothing_window_s=0.5,
                    fps=None,
                    speed_col='Speed'):
    """
      1. Smooth the distance trace to remove noise.
      2. Find 'arrival' events: local minima where the animal actually
         reached close proximity to the port (< arrival_threshold_cm).
      3. For each arrival, trace BACKWARDS in time to find the approach
         onset — defined as the last point before arrival where distance
         was at a local maximum (i.e. where the animal started moving
         toward the port without turning back more than
         max_distance_increase_cm).
      4. Flag the full onset→arrival segment as an approach bout.

    This correctly captures approaches that begin from any distance,
    ignores transient distance fluctuations, and does not produce
    spurious bouts.

    Parameters
    ----------
    arrival_threshold_cm : float
        Maximum nose-to-port distance (cm) at approach endpoint.
        The animal must actually arrive close to the port for the
        trajectory to count as an approach.
    min_approach_duration_s : float
        Minimum duration (seconds) of the approach segment.
        Filters out incidental passes where the animal briefly comes close.
    max_distance_increase_cm : float
        Maximum allowed transient increase in distance during a backward
        trace before the onset is declared. Tolerates small deviations
        (sniffing, head turning) without prematurely ending the bout.
    smoothing_window_s : float
        Gaussian smoothing window (seconds) applied to distance trace
        before onset detection. Does not affect the raw distance columns.
    fps : float or None
        Frames per second. Inferred from Time(s) if None.

    Returns
    -------
    df with columns:
        'Dist_to_airpuff_left_cm'
        'Dist_to_airpuff_right_cm'
        'Dist_to_nearest_airpuff_cm'
        'Dist_to_lickport_cm'
        'Approach_airpuff_left'
        'Approach_airpuff_right'
        'Approach_any_airpuff'
        'Approach_lickport'
    """
    from scipy.ndimage import uniform_filter1d
    from scipy.signal import argrelmin

    scale = scale_and_coords["Scale_cm_per_px"]

    # ── Infer fps ─────────────────────────────────────────────────────────────
    if fps is None:
        dt  = np.diff(df['Time(s)'].values)
        fps = 1.0 / np.median(dt[dt > 0])

    min_frames     = int(min_approach_duration_s * fps)
    smooth_frames  = max(3, int(smoothing_window_s * fps))

    # ── Raw distances (always stored) ─────────────────────────────────────────
    dist_left  = compute_distance(df, ports['airpuff_left'],  scale)
    dist_right = compute_distance(df, ports['airpuff_right'], scale)
    dist_lick  = compute_distance(df, ports['lick_port'],     scale)

    df['Dist_to_airpuff_left_cm']    = dist_left
    df['Dist_to_airpuff_right_cm']   = dist_right
    df['Dist_to_nearest_airpuff_cm'] = np.minimum(dist_left, dist_right)
    df['Dist_to_lickport_cm']        = dist_lick

    # ── Core approach detection ───────────────────────────────────────────────
    def _detect_approach_bouts(raw_dist):
        """
        Returns a binary array (1 = approaching) for a single port distance trace.
        """
        n      = len(raw_dist)
        result = np.zeros(n, dtype=int)

        # Smooth to find structure without noise
        smooth = uniform_filter1d(raw_dist.values, size=smooth_frames)

        # Step 1: find arrivals — frames where smoothed distance is a local
        # minimum AND below arrival_threshold_cm
        local_mins = argrelmin(smooth, order=min_frames)[0]
        arrivals   = local_mins[smooth[local_mins] < arrival_threshold_cm]

        if len(arrivals) == 0:
            return result

        # Step 2: for each arrival, trace backwards to find onset
        for arrival in arrivals:
            dist_at_arrival = smooth[arrival]
            max_dist_seen   = dist_at_arrival
            onset           = arrival  # will be pushed back

            for i in range(arrival - 1, -1, -1):
                d = smooth[i]

                if d > max_dist_seen + max_distance_increase_cm:
                    # Distance increased too much — approach started after this
                    onset = i + 1
                    break

                max_dist_seen = max(max_dist_seen, d)
                onset = i  # keep extending onset backwards

            # Step 3: apply minimum duration filter
            if (arrival - onset) >= min_frames:
                result[onset:arrival + 1] = 1

        return result

    # ── Apply to each port ────────────────────────────────────────────────────
    port_map = {
        'airpuff_left'  : ('Dist_to_airpuff_left_cm',  'Approach_airpuff_left'),
        'airpuff_right' : ('Dist_to_airpuff_right_cm', 'Approach_airpuff_right'),
        'lick_port'     : ('Dist_to_lickport_cm',      'Approach_lickport'),
    }

    for port_key, (dist_col, approach_col) in port_map.items():
        bouts = _detect_approach_bouts(df[dist_col])
        df[approach_col] = bouts
        n = bouts.sum()
        print(f"  [approach] {port_key}: {n} frames flagged "
              f"({n/len(df)*100:.1f}%), "
              f"{(np.diff(np.r_[0, bouts, 0]) == 1).sum()} bouts")

    df['Approach_any_airpuff'] = (
        (df['Approach_airpuff_left'] == 1) | (df['Approach_airpuff_right'] == 1)
    ).astype(int)

    return df

def detect_head_orientation(df, ports, scale_and_coords,
                             angle_threshold_deg=10.0,
                             nose_x='nose_x', nose_y='nose_y',
                             ear_left_x='left_ear_x', ear_left_y='left_ear_x',
                             ear_right_x='right_ear_x', ear_right_y='right_ear_y'):
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
                   nose_center_threshold_cm=4.0,
                   orientation_angle_threshold_deg=10.0,
                   nose_x='nose_x', nose_y='nose_y',
                   center_x='center_x', center_y='center_y',
                   ear_left_x='left_ear_x', ear_left_y='left_ear_x',
                   ear_right_x='right_ear_x', ear_right_y='right_ear_y'):
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

def categorize_airpuff_events(df, ports, scale_and_coords,
                               airpuff_col='Airpuffs',
                               pre_window_s=3.0,
                               post_window_s=3.0,
                               escape_distance_threshold_cm=5.0,
                               nose_angle_threshold_deg=45.0,
                               side_angle_threshold_deg=90.0,
                               fps=None,
                               nose_x='nose_x', nose_y='nose_y',
                               ear_left_x='ear_left_x', ear_left_y='ear_left_y',
                               ear_right_x='ear_right_x', ear_right_y='ear_right_y',
                               speed_col='Speed'):
        
    """
    Classify each airpuff event along three independent dimensions:

    DIMENSION 1 — Body part sensed
    --------------------------------
    Based on head orientation AT THE MOMENT of airpuff delivery:
      'Nose'       : head pointed toward port (angle < nose_angle_threshold_deg)
      'Head side'  : head angled toward port but not directly
                     (nose_angle_threshold_deg <= angle < side_angle_threshold_deg)
      'Flank/other': body not oriented toward port

    DIMENSION 2 — Prior assessment
    --------------------------------
    Based on mouse behaviour in the pre_window_s seconds BEFORE airpuff:
      'Assessed'    : both Head_toward AND Stretch_toward port detected
      'Looked only' : Head_toward detected but no Stretch
      'Approached'  : Approach_any_airpuff detected but no head orientation
      'Unexpected'  : none of the above — airpuff caught mouse off guard

    DIMENSION 3 — Escape behaviour
    --------------------------------
    Based on nose-to-port distance in the post_window_s seconds AFTER airpuff:
      'Active escape' : mean distance in post-window significantly exceeds
                        distance at airpuff onset (mouse moved away)
      'Stayed'        : mouse remained close to the port

    Parameters
    ----------
    df : pd.DataFrame
        Must contain Time(s), nose coordinates, and all columns produced by
        detect_head_orientation(), detect_approach(), detect_stretch().
    ports : dict
        Port coordinates with keys 'airpuff_left' and 'airpuff_right'.
    scale_and_coords : dict
        Must contain 'Scale_cm_per_px'.
    airpuff_col : str
        Column marking airpuff delivery frames (1 = delivered, 0 = not).
    pre_window_s : float
        Seconds before airpuff onset to assess prior behaviour.
    post_window_s : float
        Seconds after airpuff onset to assess escape.
    escape_distance_threshold_cm : float
        Minimum increase in nose-to-port distance (cm) between onset and
        post-window mean to classify as 'Active escape'.
    nose_angle_threshold_deg : float
        Angle cutoff for 'Nose' body part category.
    side_angle_threshold_deg : float
        Angle cutoff for 'Head side' body part category.
    fps : float or None
        Frames per second. If None, inferred from Time(s) column.

    Returns
    -------
    events_df : pd.DataFrame
        One row per airpuff event with columns:
            'Event_idx'           — index of onset frame in df
            'Time_onset_s'        — time of airpuff onset
            'Port'                — which port was closest ('left' or 'right')
            'Angle_at_onset_deg'  — head-to-port angle at onset
            'Dist_at_onset_cm'    — nose-to-port distance at onset
            'Body_part_sensed'    — Dimension 1 label
            'Prior_assessment'    — Dimension 2 label
            'Escape_behaviour'    — Dimension 3 label
            'Mean_dist_post_cm'   — mean distance in post-window (useful for QC)
            'Max_speed_post_cms'  — max speed in post-window
    """

    import warnings
    scale = scale_and_coords["Scale_cm_per_px"]

    # ── Infer fps ─────────────────────────────────────────────────────────────
    if fps is None:
        dt  = np.diff(df['Time(s)'].values)
        fps = 1.0 / np.median(dt[dt > 0])
        print(f"  [airpuff] Inferred fps: {fps:.2f} Hz")

    pre_frames  = int(pre_window_s  * fps)
    post_frames = int(post_window_s * fps)

    # ── Detect airpuff onset frames (0→1 transitions) ─────────────────────────
    airpuff_vals = df[airpuff_col].fillna(0).values.astype(int)
    onsets = np.where(np.diff(np.r_[0, airpuff_vals]) == 1)[0]
    print(f"  [airpuff] {len(onsets)} airpuff events detected.")

    if len(onsets) == 0:
        warnings.warn("No airpuff events found. Check airpuff_col and its encoding.")
        return pd.DataFrame()

    # ── Precompute distances and angles for both ports ────────────────────────
    dist_left  = compute_distance(df, ports['airpuff_left'],  scale,
                                  nose_x=nose_x, nose_y=nose_y).values
    dist_right = compute_distance(df, ports['airpuff_right'], scale,
                                  nose_x=nose_x, nose_y=nose_y).values

    # Head orientation angles (may already be in df from detect_head_orientation)
    def _get_or_compute_angle(port_key, angle_col):
        if angle_col in df.columns:
            return df[angle_col].values
        # Recompute on the fly if column missing
        ear_mid_x  = (df[ear_left_x]  + df[ear_right_x])  / 2
        ear_mid_y  = (df[ear_left_y]  + df[ear_right_y]) / 2
        head_vec_x = df[nose_x] - ear_mid_x
        head_vec_y = df[nose_y] - ear_mid_y
        port_vec_x = ports[port_key]['x'] - df[nose_x]
        port_vec_y = ports[port_key]['y'] - df[nose_y]
        dot = head_vec_x * port_vec_x + head_vec_y * port_vec_y
        mag = (np.sqrt(head_vec_x**2 + head_vec_y**2) *
               np.sqrt(port_vec_x**2 + port_vec_y**2) + 1e-8)
        return np.degrees(np.arccos(np.clip(dot / mag, -1.0, 1.0)))

    angles_left  = _get_or_compute_angle('airpuff_left',  'Angle_to_airpuff_left')
    angles_right = _get_or_compute_angle('airpuff_right', 'Angle_to_airpuff_right')

    # ── Optional pre-computed behaviour columns ────────────────────────────────
    def _get_col(name):
        return df[name].values if name in df.columns else np.zeros(len(df), dtype=int)

    head_toward  = _get_col('Head_toward_any_airpuff')
    stretch      = _get_col('Stretch_toward_any_airpuff')
    approach     = _get_col('Approach_any_airpuff')
    speed        = df[speed_col].values if speed_col in df.columns else np.zeros(len(df))

    # ── Classify each event ───────────────────────────────────────────────────
    records = []

    for onset in onsets:
        # Pre- and post-window frame ranges (clamped to df bounds)
        pre_start  = max(0, onset - pre_frames)
        post_end   = min(len(df) - 1, onset + post_frames)

        # Which port is closest at onset?
        if dist_left[onset] <= dist_right[onset]:
            port_label  = 'left'
            dist_onset  = dist_left[onset]
            angle_onset = angles_left[onset]
            dist_post   = dist_left[pre_start:post_end]      # reuse for post
            dist_post   = dist_left[onset:post_end]
        else:
            port_label  = 'right'
            dist_onset  = dist_right[onset]
            angle_onset = angles_right[onset]
            dist_post   = dist_right[onset:post_end]

        # ── Dimension 1: Body part sensed ─────────────────────────────────────
        if angle_onset < nose_angle_threshold_deg:
            body_part = 'Nose'
        elif angle_onset < side_angle_threshold_deg:
            body_part = 'Head side'
        else:
            body_part = 'Flank/other'

        # ── Dimension 2: Prior assessment ─────────────────────────────────────
        pre_head    = head_toward[pre_start:onset].sum()
        pre_stretch = stretch[pre_start:onset].sum()
        pre_approach= approach[pre_start:onset].sum()

        if pre_head > 0 and pre_stretch > 0:
            prior = 'Assessed'
        elif pre_head > 0:
            prior = 'Looked only'
        elif pre_approach > 0:
            prior = 'Approached'
        else:
            prior = 'Unexpected'

        # ── Dimension 3: Escape behaviour ─────────────────────────────────────
        mean_dist_post = float(np.mean(dist_post)) if len(dist_post) > 0 else np.nan
        max_speed_post = float(np.max(speed[onset:post_end])) if post_end > onset else np.nan

        if mean_dist_post - dist_onset > escape_distance_threshold_cm:
            escape = 'Active escape'
        else:
            escape = 'Stayed'

        records.append({
            'Event_idx'          : onset,
            'Time_onset_s'       : df['Time(s)'].iloc[onset],
            'Port'               : port_label,
            'Angle_at_onset_deg' : round(angle_onset, 2),
            'Dist_at_onset_cm'   : round(dist_onset, 2),
            'Body_part_sensed'   : body_part,
            'Prior_assessment'   : prior,
            'Escape_behaviour'   : escape,
            'Mean_dist_post_cm'  : round(mean_dist_post, 2),
            'Max_speed_post_cms' : round(max_speed_post, 2),
        })

    events_df = pd.DataFrame(records)

    # ── Summary printout ──────────────────────────────────────────────────────
    print(f"\n  ── Airpuff event summary ({len(events_df)} events) ──")
    for dim, col in [('Body part sensed', 'Body_part_sensed'),
                     ('Prior assessment', 'Prior_assessment'),
                     ('Escape behaviour', 'Escape_behaviour')]:
        print(f"\n  {dim}:")
        for label, count in events_df[col].value_counts().items():
            print(f"    {label:<20} {count:>3} ({count/len(events_df)*100:.1f}%)")

    return events_df