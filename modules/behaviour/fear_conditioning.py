import pandas as pd
import re
import numpy as np

def parse_protocol_sheet(path, sheet_name):
    """
    Parse protocol excel sheet and extract intervals for:
      - CS+ (SND column 'on(F1)' ... '!on')
      - CS- (SND column 'on(F3)' ... '!on')
      - Shock (LED2(1,3) column 'ON' ... '!on')
      - LED3 (LED3(1,4) column 'on' ... '!on') for protocol-wide on/off
    
    Returns dict with lists of (start_s, end_s) in seconds (relative to protocol start).
    """
    df = pd.read_excel(path, sheet_name=sheet_name, dtype=str).fillna("")
    # Normalize columns (strip whitespaces)
    df.columns = [str(c).strip() for c in df.columns]

    time_col = 'T1'
    snd_col = 'SND(6,2)'
    shk_col = 'LED2(1,3)'
    led3_col = 'LED3(1,4)'
    
    # Convert ms → seconds
    df[time_col] = df[time_col] / 1000.0

    # Generic extractor that looks for 'on' and '!on' tokens (case-insensitive),
    # but with optional extras like on(F1), on(F3)
    def extract_signal_intervals(df, col, on_match_fn):
        """
        on_match_fn(cell_text_lower) -> returns:
            - 'on' if row starts or indicates ON,
            - '!on' if row indicates OFF,
            - '' otherwise
        """
        intervals = []
        current_start = None
        running_time = 0.0
        for _, row in df.iterrows():
            duration = float(row[time_col])
            cell = str(row.get(col, "")).strip()

            flag = on_match_fn(cell.lower())

            if flag == 'on':
                if current_start is None:
                    current_start = running_time
            elif flag == '!on':
                # OFF at beginning of this row => interval closes at running_time
                if current_start is not None:
                    intervals.append((current_start, running_time))
                    current_start = None
            # else no change

            running_time += duration

        # close if file ends while ON
        if current_start is not None:
            intervals.append((current_start, running_time))

        return intervals

    # Match functions
    def snd_match(cell_lower):
        # return 'on' for on(f1) or on(f3), distinguishing later outside
        # but we only want to detect on(f1) and on(f3) here; for general 'on' in SND ignore others
        if re.search(r'on\s*\(\s*f1\s*\)', cell_lower):
            return 'on_f1'
        if re.search(r'on\s*\(\s*f3\s*\)', cell_lower):
            return 'on_f3'
        if '!on' in cell_lower:
            return '!on'
        return ''

    # We'll run through SND column row-by-row and build intervals for f1 and f3
    cs_plus = []
    cs_minus = []
    running_time = 0.0
    # we assume rows in SND can contain either "on(F1)", "!on", or "on(F3)", etc.
    current_f1 = None
    current_f3 = None
    for _, row in df.iterrows():
        duration = float(row[time_col])
        cell = str(row.get(snd_col, "")).strip().lower()

        # check for explicit on(F1) or on(F3)
        if re.search(r'on\s*\(\s*f1\s*\)', cell):
            if current_f1 is None:
                current_f1 = running_time
        if re.search(r'on\s*\(\s*f3\s*\)', cell):
            if current_f3 is None:
                current_f3 = running_time

        # check for !on: closes whichever sound(s) are currently open
        if '!on' in cell:
            if current_f1 is not None:
                cs_plus.append((current_f1, running_time))
                current_f1 = None
            if current_f3 is not None:
                cs_minus.append((current_f3, running_time))
                current_f3 = None

        running_time += duration

    # Close any remaining open ones at end of sheet
    total_dur = df[time_col].sum()
    if current_f1 is not None:
        cs_plus.append((current_f1, total_dur))
    if current_f3 is not None:
        cs_minus.append((current_f3, total_dur))

    # Shock extraction using same on/off logic as earlier extractor, but shock appears as 'ON' and '!on'
    def shk_match(cell_lower):
        if 'on' == cell_lower.strip().lower() or 'on' in cell_lower:
            # treat any 'on' as shock on (you can tighten logic if needed)
            return 'on'
        if '!on' in cell_lower:
            return '!on'
        return ''

    # Use the generic extractor for shock and led3
    # For shock we treat exact cell match having 'on' text; some rows use 'ON' or 'on(F1)' for other columns so be careful
    def generic_on_fn(cell_lower):
        if 'on' in cell_lower and '!on' not in cell_lower:
            return 'on'
        if '!on' in cell_lower:
            return '!on'
        return ''

    shock = extract_signal_intervals(df, shk_col, generic_on_fn)
    led3 = extract_signal_intervals(df, led3_col, generic_on_fn)

    # Done
    return {
        "CS+": cs_plus,
        "CS-": cs_minus,
        "Shock": shock,
        "LED3": led3
    }
