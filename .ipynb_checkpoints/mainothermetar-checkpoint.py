# download_metar_iem.py

import argparse
import os
import requests
from datetime import datetime

def download_iem_metar(station, start_date, end_date, out_dir='output'):
    os.makedirs(out_dir, exist_ok=True)
    
    # Iowa State University METAR request URL
    base_url = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
    
    # Format output filename
    filename = f"{station}_{start_date}_to_{end_date}.csv"
    filepath = os.path.join(out_dir, filename)

    # Split date components
    y1, m1, d1 = start_date.split("-")
    y2, m2, d2 = end_date.split("-")

    params = {
        "station": station,
        "data": ["tmpf", "dwpf", "relh", "drct", "sknt", "p01i", "alti", "vsby", "gust"],
        "year1": y1,
        "month1": m1,
        "day1": d1,
        "year2": y2,
        "month2": m2,
        "day2": d2,
        "tz": "Etc/UTC",
        "format": "comma",
        "latlon": "yes",
        "missing": "M",
        "trace": "T",
        "direct": "yes",
        "report_type": "3",  # 1: METAR, 2: SPECI, 3: both
    }

    print(f"🔍 Requesting data for station {station} from {start_date} to {end_date}...")

    response = requests.get(base_url, params=params)

    if response.status_code == 200 and "station" in response.text.lower():
        with open(filepath, "w") as f:
            f.write(response.text)
        print(f"✅ Data saved to {filepath}")
        return filepath
    else:
        print("❌ Failed to retrieve data.")
        print("Status code:", response.status_code)
        print("Response preview:", response.text[:300])
        return None

# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download METAR data from IEM archive")
    parser.add_argument("--station", required=True, help="4-letter station code (e.g., KDEN)")
    parser.add_argument("--start", required=True, help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end", required=True, help="End date (YYYY-MM-DD)")
    args = parser.parse_args()

    download_iem_metar(args.station, args.start, args.end)

