
# open ecmwf data available for download since January 21 2022 00:00 ('2022-01-21 00:00')

import argparse
from src.Ejemplo_8 import plot_nwp


def main (date_nwp,extent):
       plot_nwp(model='ecmwf',
         product='oper',
         #date_nwp='2024-01-23 12:00',
         date_nwp=date_nwp,       
         area='custom',
         #extent=[-140.0, -70.00, 70.00, 70.00],
         extent=extent,
         fxx='date',
         var='q', scale=10**3, offset=0,
         level_hpa='850',
         level_min=0, level_max=15, level_int=2,
         plot_type='c_unfilled', contour_color='black',
         apply_cmap=True, cmap='rainbow_r',
         view_clabel=True, fontsize=6,
         linewidth=0.8, linestyle='solid',
         land_ocean=True, land_color='lightgray', ocean_color='white',
         alpha=1.0, title='new',
         figsize=[15,15])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Humidity using ECMWF data")
    parser.add_argument("--date_nwp", type=str, required=True, help="NWP date-time, e.g. '2023-11-01 12:00'")
    parser.add_argument("--extent", type=float, nargs=4, required=True, metavar=('WEST', 'SOUTH', 'EAST', 'NORTH'), help="Map extent: west south east north")
    args = parser.parse_args()
    main(args.date_nwp, args.extent)

