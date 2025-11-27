import pandas as pd

def parse_protocol_sheet(path, sheet_name):
    """
    Parse protocol excel sheet and extract intervals for:
      - CS+ (CS+ column 'on' ... '!on')
      - CS- (CS- column 'on' ... '!on')
      - Shock (LED2(1,3) column 'ON' ... '!on')
      - LED3 (LED3(1,4) column 'on' ... '!on') for protocol-wide on/off
    
    Returns dict with lists of (start_s, end_s) in seconds (relative to protocol start).
    """
    df = pd.read_excel(path, sheet_name=sheet_name, dtype=str).fillna("")
    # Normalize columns (strip whitespaces)
    df.columns = [str(c).strip() for c in df.columns]

    time_col = 'T1'
    csplus_col = 'CS+'
    csminus_col = 'CS-'
    shk_col = 'LED2(1,3)'
    led3_col = 'LED3(1,4)'
    
    # Convert ms → seconds
    df[time_col] = df[time_col] / 1000.0

    # ───────────────────────────────────────────────
    # Helper: extract ON/OFF intervals for a column
    # ───────────────────────────────────────────────
    def extract_intervals(df, signal_col):
        intervals = []
        current_start = None
        current_time = 0.0  # running protocol time
        
        for idx, row in df.iterrows():
            duration = row[time_col]
            cell = str(row.get(signal_col, "")).lower()

            if 'on' == cell:          # ON starts now
                if current_start is None:
                    current_start = current_time
                    
            elif '!on' in cell:        # OFF at start of this row
                if current_start is not None:
                    intervals.append((current_start, current_time))
                    current_start = None
            
            current_time += duration

        # Close if sheet ends while ON:
        if current_start is not None:
            intervals.append((current_start, current_time))
        
        return intervals

    # Extract LED3 → defines protocol start/stop
    led3 = extract_intervals(df, led3_col)

    # Extract CS+ / CS– from columns
    cs_plus = extract_intervals(df, csplus_col)
    cs_minus = extract_intervals(df, csminus_col)
    
    # Extract shock intervals
    shock = extract_intervals(df, shk_col)

    return {
        "CS+": cs_plus,
        "CS-": cs_minus,
        "Shock": shock,
        "LED3": led3
    }
