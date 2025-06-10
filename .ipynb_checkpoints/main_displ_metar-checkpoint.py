#UNIDATA's THREDDS TDS is designed primarily for real-time and recent data access, not long-term archiving. Historical METAR archives are not retained long-term on this server- ~7 days (rolling)
# main_displ_metar.py

# 
import argparse
from src.Ejemplo_13 import plot_metar  # Make sure this file and function exist

#def main():
def main (date_metar,extent):
    plot_metar(date_metar=date_metar,
        area='custom',
        #extent=[-140.0, -70.0, 70.0, 70.0],
        extent=extent,
        radius=3, land_ocean=True, land_color = 'beige', ocean_color = 'skyblue',
        figsize=[15, 15])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MSLP using ECMWF data")
    parser.add_argument("--date_metar", type=str, required=True, help="Metar date-time, e.g. '2023-11-01 12:00'")
    parser.add_argument("--extent", type=float, nargs=4, required=True, metavar=('WEST', 'SOUTH', 'EAST', 'NORTH'), help="Map extent: west south east north")
    args = parser.parse_args()
    main(args.date_metar, args.extent)

