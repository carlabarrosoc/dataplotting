DATE_NWP="2025-05-26 00:00"
#WEST SOUTH EAST NORTH
EXTENT="-20.0 -40.0 60.0 40.0"

python main_displ_geopot.py --date_nwp "$DATE_NWP" --extent $EXTENT
python main_displ_mslp.py --date_nwp "$DATE_NWP" --extent $EXTENT
python main_displ_humid.py --date_nwp "$DATE_NWP" --extent $EXTENT
python main_displ_windbarbs.py --date_nwp "$DATE_NWP" --extent $EXTENT
python main_displ_windspeed.py --date_nwp "$DATE_NWP" --extent $EXTENT
python main_displ_metar.py --date_metar "$DATE_NWP" --extent $EXTENT

