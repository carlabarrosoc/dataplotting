
# open ecmwf data available for download since January 21 2022 00:00 ('2022-01-21 00:00')
# gfs data available for download since February 26 2021 00:00 ('2021-02-26 00:00')
# gfs products: 'pgrb2.0p25', 'pgrb2.0p50', 'pgrb2.1p00'

# geopotential height (500 hPa [4800~6000] and 300 hPa [8300~9700])

import argparse
from Ejemplo_8 import plot_nwp


def main (date_nwp,extent):
       plot_nwp(model='ecmwf',
         product='oper',
         #date_nwp='2024-01-23 12:00',
         date_nwp=date_nwp,       
         area='custom',
         #extent=[-140.0, -70.00, 70.00, 70.00],
         extent=extent,
         fxx='date',
         var='gh', scale=1, offset=0,
         level_hpa='500',
         level_min=4800, level_max=6000, level_int=60,
         plot_type='c_unfilled', contour_color='black',
         apply_cmap=True, cmap='turbo',
         view_clabel=True, fontsize=6,
         linewidth=0.8, linestyle='solid',
         land_ocean=True, land_color='lightgray', ocean_color='white',
         alpha=1.0, title='new',
         figsize=[15,15])

#if __name__ == "__main__":
#    main()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MSLP using ECMWF data")
    parser.add_argument("--date_nwp", type=str, required=True, help="NWP date-time, e.g. '2023-11-01 12:00'")
    parser.add_argument("--extent", type=float, nargs=4, required=True, metavar=('WEST', 'SOUTH', 'EAST', 'NORTH'), help="Map extent: west south east north")
    args = parser.parse_args()
    main(args.date_nwp, args.extent)

