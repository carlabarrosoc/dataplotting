
# Mean Sea Level Pressure (MSLP) plotting using ECMWF data
import argparse
from src.Ejemplo_8 import plot_nwp  # Make sure this file and function exist

#def main():
def main (date_nwp,extent):
    plot_nwp(model='ecmwf',
        product='oper',
        #date_nwp='2024-01-01 12:00',
        date_nwp=date_nwp,
        area='custom',
        #extent=[-140.0, -70.0, 70.0, 70.0],
        extent=extent,
        fxx='date',  # forecast hour, if needed
        var='windspeed',scale=1,offset=0,
        level_hpa='500',
        level_min=-0,level_max=70,level_int=5,
        plot_type='c_filled',contour_color='gray',
        apply_cmap=True, cmap='GnBu',
        view_clabel=True, fontsize=6,
        linewidth=0.5, linestyle='solid',
        land_ocean=False, land_color='lightgray', ocean_color='white',
        alpha=1.0, title='new',
        figsize=[15, 15])
  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MSLP using ECMWF data")
    parser.add_argument("--date_nwp", type=str, required=True, help="NWP date-time, e.g. '2023-11-01 12:00'")
    parser.add_argument("--extent", type=float, nargs=4, required=True, metavar=('WEST', 'SOUTH', 'EAST', 'NORTH'), help="Map extent: west south east north")
    args = parser.parse_args()
    main(args.date_nwp, args.extent)