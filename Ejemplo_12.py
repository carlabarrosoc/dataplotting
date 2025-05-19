"""
========================================================================
METEOROLOGICAL DATA VISUALIZATION - EJEMPLO_12.PY
========================================================================

PURPOSE:
    This script visualizes numerical weather prediction (NWP) data with
    an enhanced focus on identifying and marking pressure systems.
    It generates maps with high/low pressure markers and contours.

CREATED BY:
    Carla Barroso, EUMETSAT
    Last Updated: 2025-05-19

    A portion of this work used code generously provided by Brian Blaylock's 
    Carpenter Workshop python package (https://github.com/blaylockbk/Carpenter_Workshop)

USAGE:
    Run this script directly with Python:
    $ python Ejemplo_12.py
    
    The script uses default parameters but these can be modified in the 
    script execution section at the bottom of the file.

REQUIRED LIBRARIES:
    - cartopy: For geographical plotting
    - matplotlib: For visualization
    - herbie: For accessing weather model data
    - metpy: For meteorological calculations
    - netCDF4: For reading data files
"""

#===========================================================================
# IMPORT LIBRARIES
#===========================================================================

import os
import glob
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time
from datetime import timedelta, date, datetime
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from herbie import Herbie
from metpy.io import metar
from metpy.plots import current_weather, sky_cover, StationPlot
from metpy.calc import reduce_point_density
from netCDF4 import Dataset

#===========================================================================
# MAP PROJECTION SETUP
#===========================================================================

pc = ccrs.PlateCarree()
pc._threshold = 0.01  # https://github.com/SciTools/cartopy/issues/8

#===========================================================================
# MAP DRAWING CLASS
#===========================================================================

class EasyMap:
    """
    Build a matplotlib/cartopy axes with common map elements.

    This class does about 95% of my Cartopy needs.

    Feature elements from https://www.naturalearthdata.com/features/

    TODO: Rename to CommonFeatures or EasyMap
    """

    def __init__(
        self,
        scale="110m",
        ax=None,
        crs=pc,
        *,
        figsize=None,
        fignum=None,
        dpi=None,
        dark=False,
        verbose=False,
        add_coastlines=True,
        facecolor=None,
        coastlines_kw={},
        **kwargs,
    ):
        """
        Initialize a Cartopy axes

        Add coastlines with this method. Use other methods to add
        other common features to the Cartopy axes.

        .. tip:: ``ax=None`` is a great way to initialize a new Cartopy axes.

        Methods
        -------
        .adjust_extent
        .center_extent
        .copy_extent

        Parameters
        ----------
        scale : {'10m', '50m' 110m'}
            The cartopy feature's level of detail.
            .. note::  The ``'10m'`` scale for OCEAN and LAND takes a *long* time.
        ax : plot axes
            The axis to add the feature to.
            If None, it will create a new cartopy axes with ``crs``.
        crs : cartopy.crs
            Coordinate reference system (aka "projection") to create new map
            if no cartopy axes is given. Default is ccrs.PlateCarree.
        dark : bool
            If True, use alternative "dark theme" colors for land and water.

            .. figure:: _static/BB_maps/common_features-1.png
            .. figure:: _static/BB_maps/common_features-2.png

        add_coastlines : bool
            For convince, the coastlines are added to the axes by default.
            This can be turned off and instead controlled with the COASTLINES
            method.
        coastlines_kw : dict
            kwargs for the default COASTLINES method.

        figsize : tuple or float
            Set the figure size.
            If single number given, then will make a square figure.
        fignum : int
            If given, create a new figure and supblot for the given crs.
            (This might be handy in a loop when you want to create maps on
            several figures instead of plotting on the same figure.)
        dpi : int
            Set the figure dpi

        Examples
        --------
        https://github.com/blaylockbk/Carpenter_Workshop/blob/main/notebooks/demo_cartopy_tools.ipynb

        >>> feat = EasyMap()
        >>> feat.OCEAN().STATES()
        >>> ax = feat.ax

        Alternatively,

        >>> ax = EasyMap().ax
        >>> feat = ax.EasyMap
        >>> feat.OCEAN().STATES()
        """
        self.scale = scale
        self.ax = ax
        self.crs = crs
        self.figsize = figsize
        self.fignum = fignum
        self.dpi = dpi
        self.dark = dark
        self.verbose = verbose
        self.kwargs = kwargs

        self.ax = check_cartopy_axes(
            ax=self.ax, crs=self.crs, fignum=self.fignum, verbose=self.verbose
        )

        # In a round-about way, you can get this EasyMap object from the axes
        # >>> ax = EasyMap().ax
        # >>> ax.EasyMap.STATES()
        self.ax.EasyMap = self

        self.kwargs.setdefault("linewidth", 0.75)

        # NOTE: I don't use the 'setdefault' method here because it doesn't
        # work as expect when switching between dark and normal themes.
        # The defaults would be set the first time the function is called,
        # but the next time it is called and `dark=True` the defaults do not
        # reset. I don't know why this is the behavior.
        if self.dark:
            self.land = "#060613"  # dark (default)
            self.land1 = "#3f3f3f"  # lighter (alternative)
            self.water = "#0f2b38"
            # https://github.com/SciTools/cartopy/issues/880
            self.ax.set_facecolor(self.land)  # requires cartopy >= 0.18
            self.kwargs = {**{"edgecolor": ".5"}, **self.kwargs}
        else:
            self.land = "#efefdb"  # tan (default)
            self.land1 = "#dbdbdb"  # grey (alternative)
            self.water = "#97b6e1"
            self.kwargs = {**{"edgecolor": ".15"}, **self.kwargs}

        if facecolor:
            # Instead of applying both LAND and OCEAN,
            # it may be faster to just set the facecolor of land
            # and then only apply the OCEAN method.
            if facecolor.lower() == "land":
                self.ax.set_facecolor(self.land)
            elif facecolor.lower() == "land1":
                self.ax.set_facecolor(self.land1)
            elif facecolor.lower() == "water":
                self.ax.set_facecolor(self.water)
            else:
                self.ax.set_facecolor(facecolor)

        if add_coastlines:
            # Default map will automatically add COASTLINES
            self.COASTLINES(**coastlines_kw)

        if figsize is not None:
            if hasattr(figsize, "__len__"):
                plt.gcf().set_figwidth(self.figsize[0])
                plt.gcf().set_figheight(self.figsize[1])
            else:
                plt.gcf().set_figwidth(self.figsize)
                plt.gcf().set_figheight(self.figsize)
        if dpi is not None:
            plt.gcf().set_dpi(self.dpi)

        # Add my custom methods
        self.ax.__class__.adjust_extent = _adjust_extent
        self.ax.__class__.center_extent = _center_extent
        self.ax.__class__.copy_extent = _copy_extent

    # Feature Elements
    def COASTLINES(self, **kwargs):
        kwargs.setdefault("zorder", 100)
        kwargs.setdefault("facecolor", "none")
        kwargs = {**self.kwargs, **kwargs}
        self.ax.add_feature(feature.COASTLINE.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 COASTLINES:", kwargs)
        return self

    def BORDERS(self, **kwargs):
        """Borders between countries. *Excludes coastlines*"""
        kwargs.setdefault("linewidth", 0.5)
        kwargs = {**self.kwargs, **kwargs}
        self.ax.add_feature(feature.BORDERS.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 BORDERS:", kwargs)
        return self

    def STATES(self, **kwargs):
        """State and Province borders. *Includes coastlines*

        Note: If scale="110m", only the US States are drawn.
              If scale="50m", then more country states/provinces are drawn.
              If scale="10m", then even *more* countries drawn.
        """
        kwargs.setdefault("alpha", 0.15)

        kwargs = {**self.kwargs, **kwargs}
        self.ax.add_feature(feature.STATES.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 STATES:", kwargs)
        return self

    def STATES2(self, **kwargs):
        """States and Provinces (US, Canada, Australia, Brazil, China, Inda, etc.)

        Alternative source for data than provided by STATES.
        """
        kwargs.setdefault("alpha", 0.15)

        kwargs = {**self.kwargs, **kwargs}
        states_provinces = feature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lines",
            scale="50m",
            facecolor="none",
        )
        self.ax.add_feature(states_provinces, **kwargs)

        if self.verbose == "debug":
            print("🐛 STATES2:", kwargs)
        return self

    def COUNTIES(self, counties_scale="20m", **kwargs):
        """US counties. *Includes coastslines*"""
        _counties_scale = {"20m", "5m", "500k"}
        assert (
            counties_scale in _counties_scale
        ), f"counties_scale must be {_counties_scale}"
        kwargs.setdefault("linewidth", 0.33)
        kwargs.setdefault("alpha", 0.15)
        kwargs = {**self.kwargs, **kwargs}
        self.ax.add_feature(USCOUNTIES.with_scale(counties_scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 COUNTIES:", kwargs)
        return self

    def OCEAN(self, **kwargs):
        """Color-filled ocean area"""
        kwargs.setdefault("edgecolor", "none")
        kwargs = {**self.kwargs, **kwargs}

        if self.dark:
            kwargs = {**{"facecolor": self.water}, **kwargs}

        self.ax.add_feature(feature.OCEAN.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 OCEAN:", kwargs)
        return self

    def LAND(self, **kwargs):
        """Color-filled land area"""
        kwargs.setdefault("edgecolor", "none")
        kwargs.setdefault("linewidth", 0)
        kwargs = {**self.kwargs, **kwargs}

        if self.dark:
            kwargs = {**{"facecolor": self.land}, **kwargs}

        self.ax.add_feature(feature.LAND.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 LAND:", kwargs)
        return self

    def RIVERS(self, **kwargs):
        """Rivers"""
        kwargs.setdefault("linewidth", 0.3)
        kwargs = {**self.kwargs, **kwargs}

        if self.dark:
            kwargs = {**{"color": self.water}, **kwargs}
        else:
            kwargs = {**{"color": self.water}, **kwargs}

        self.ax.add_feature(feature.RIVERS.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 RIVERS:", kwargs)
        return self

    def LAKES(self, **kwargs):
        """Color-filled lake area"""
        kwargs.setdefault("linewidth", 0)
        kwargs = {**self.kwargs, **kwargs}

        if self.dark:
            kwargs = {**{"facecolor": self.water}, **kwargs}
            kwargs = {**{"edgecolor": self.water}, **kwargs}
        else:
            kwargs = {**{"facecolor": feature.COLORS["water"]}, **kwargs}
            kwargs = {**{"edgecolor": feature.COLORS["water"]}, **kwargs}

        self.ax.add_feature(feature.LAKES.with_scale(self.scale), **kwargs)
        if self.verbose == "debug":
            print("🐛 LAKES:", kwargs)
        return self

    def TERRAIN(
        self,
        coarsen=30,
        *,
        top="ice",
        kind="pcolormesh",
        extent=None,
        **kwargs,
    ):
        """
        Plot terrain data from ETOPO1 dataset.

        Parameters
        ----------
        coarsen : int
            ETOPO1 data is a 1-minute arc dataset. This is huge.
            For global plots, you don't need this resolution, and can
            be happy with a 30-minute arc resolution (default).
        top : {"ice", "bedrock"}
            Top of the elevation model. "ice" is top of ice sheets in
            Greenland and Antarctica and "bedrock" is elevation of
            of ground under the ice.
        kind : {"contourf", "pcolormesh"}
            Plot data as a contour plot or pcolormesh
        extent :
            Trim the huge dataset to a specific region. (Variable cases).
            - by hemisphere {"NE", "SE", "NW", "SW"}
            - by region {"CONUS"}
            - by extent (len==4 tuple/list), e.g. `[-130, -100, 20, 50]`
            - by xarray.Dataset (must have coordinates 'lat' and 'lon')
              TODO: Currently does not allow domains that cross -180 lon.
        """
        da = get_ETOPO1(top=top, coarsen=coarsen)

        if extent:
            if isinstance(extent, (list, tuple)):
                assert (
                    len(extent) == 4
                ), "extent tuple must be len 4 (minLon, maxLon, minLat, maxLat)"
            elif isinstance(extent, str):
                assert (
                    extent in _extents
                ), f"extent string must be one of {_extents.keys()}"
                extent = _extents[extent]
            elif hasattr(extent, "coords"):
                # Get extent from lat/lon bounds in xarray DataSet
                extent = extent.rename({"latitude": "lat", "longitude": "lon"})
                extent["lon"] = _to_180(extent["lon"])
                extent = (
                    extent.lon.min().item(),
                    extent.lon.max().item(),
                    extent.lat.min().item(),
                    extent.lat.max().item(),
                )

            da = da.where(
                (da.lon >= extent[0])
                & (da.lon <= extent[1])
                & (da.lat >= extent[2])
                & (da.lat <= extent[3])
            )

        # Get "land" points (elevation is 0 and above, crude estimation)
        da = da.where(da >= 0)

        kwargs.setdefault("zorder", 0)
        kwargs.setdefault("cmap", "YlOrBr")
        kwargs.setdefault("levels", range(0, 8000, 500))
        kwargs.setdefault("vmin", 0)
        kwargs.setdefault("vmax", 8000)

        if kind == "contourf":
            _ = kwargs.pop("vmax")
            _ = kwargs.pop("vmin")
            self.ax.contourf(da.lon, da.lat, da, transform=pc, **kwargs)
        elif kind == "pcolormesh":
            _ = kwargs.pop("levels")
            self.ax.pcolormesh(da.lon, da.lat, da, transform=pc, **kwargs)

        return self

    def BATHYMETRY(
        self,
        coarsen=30,
        *,
        top="ice",
        kind="pcolormesh",
        extent=None,
        **kwargs,
    ):
        """
        Plot bathymetry data from ETOPO1 dataset.

        Parameters
        ----------
        coarsen : int
            ETOPO1 data is a 1-minute arc dataset. This is huge.
            For global plots, you don't need this resolution, and can
            be happy with a 30-minute arc resolution (default).
        top : {"ice", "bedrock"}
            Top of the elevation model. "ice" is top of ice sheets in
            Greenland and Antarctica and "bedrock" is elevation of
            of ground under the ice.
        kind : {"contourf", "pcolormesh"}
            Plot data as a contour plot or pcolormesh
        extent :
            Trim the huge dataset to a specific region. (Variable cases).
            - by hemisphere {"NE", "SE", "NW", "SW"}
            - by region {"CONUS"}
            - by extent (len==4 tuple/list), e.g. `[-130, -100, 20, 50]`
            - by xarray.Dataset (must have coordinates 'lat' and 'lon')
              TODO: Currently does not allow domains that cross -180 lon.
        """
        da = get_ETOPO1(top=top, coarsen=coarsen)

        if extent:
            if isinstance(extent, (list, tuple)):
                assert (
                    len(extent) == 4
                ), "extent tuple must be len 4 (minLon, maxLon, minLat, maxLat)"
            elif isinstance(extent, str):
                assert (
                    extent in _extents
                ), f"extent string must be one of {_extents.keys()}"
                extent = _extents[extent]
            elif hasattr(extent, "coords"):
                # Get extent from lat/lon bounds in xarray DataSet
                extent = extent.rename({"latitude": "lat", "longitude": "lon"})
                extent["lon"] = _to_180(extent["lon"])
                extent = (
                    extent.lon.min().item(),
                    extent.lon.max().item(),
                    extent.lat.min().item(),
                    extent.lat.max().item(),
                )

            da = da.where(
                (da.lon >= extent[0])
                & (da.lon <= extent[1])
                & (da.lat >= extent[2])
                & (da.lat <= extent[3])
            )

        # Get "water" points (elevation is 0 and above, crude estimation)
        da = da.where(da <= 0)

        kwargs.setdefault("zorder", 0)
        kwargs.setdefault("cmap", "Blues_r")
        kwargs.setdefault("levels", range(-10000, 1, 500))
        kwargs.setdefault("vmax", 0)
        kwargs.setdefault("vmin", -10000)

        if kind == "contourf":
            _ = kwargs.pop("vmax")
            _ = kwargs.pop("vmin")
            self.ax.contourf(da.lon, da.lat, da, transform=pc, **kwargs)
        elif kind == "pcolormesh":
            _ = kwargs.pop("levels")
            self.ax.pcolormesh(da.lon, da.lat, da, transform=pc, **kwargs)

        return self

    def PLAYAS(self, **kwargs):
        """Color-filled playa area"""
        kwargs.setdefault("linewidth", 0)
        kwargs = {**self.kwargs, **kwargs}

        if self.dark:
            kwargs = {**{"facecolor": "#4D311A73"}, **kwargs}
            kwargs = {**{"edgecolor": "none"}, **kwargs}
        else:
            kwargs = {**{"facecolor": "#FDA65473"}, **kwargs}
            kwargs = {**{"edgecolor": "none"}, **kwargs}

        playa = feature.NaturalEarthFeature("physical", "playas", "10m")
        self.ax.add_feature(playa, **kwargs)
        if self.verbose == "debug":
            print("🐛 PLAYAS:", kwargs)
        return self

    def TIMEZONE(self, **kwargs):
        """Timezone boundaries"""
        kwargs.setdefault("linewidth", 0.2)
        kwargs.setdefault("facecolor", "none")
        kwargs.setdefault("linestyle", ":")
        kwargs = {**self.kwargs, **kwargs}
        tz = feature.NaturalEarthFeature("cultural", "time_zones", "10m")
        self.ax.add_feature(tz, **kwargs)
        if self.verbose == "debug":
            print("🐛 TIMEZONE:", kwargs)
        return self

    def ROADS(self, road_types=None, **kwargs):
        """
        Major roads

        Parameters
        ----------
        road_types : None, str, list
            Filter the types of roads you want. The road type may be a single
            string or a list of road types.
            e.g. ['Major Highway', 'Secondary Highway']

        Of course, the shapefile has many other road classifiers for each
        road, like "level" (Federal, State, Interstate), road "name",
        "length_km", etc. Filters for each of these could be added if I
        need them later.

        """

        kwargs.setdefault("edgecolor", "#b30000")
        kwargs.setdefault("facecolor", "none")
        kwargs.setdefault("linewidth", 0.2)

        kwargs = {**self.kwargs, **kwargs}

        if road_types is None:
            # Plot all roadways
            roads = feature.NaturalEarthFeature("cultural", "roads", "10m", **kwargs)
            self.ax.add_feature(roads)
        else:
            # Specify the type of road to include in plot
            if isinstance(road_types, str):
                road_types = [road_types]
            shpfilename = shapereader.natural_earth("10m", "cultural", "roads")
            df = geopandas.read_file(shpfilename)
            _types = df["type"].unique()
            assert np.all(
                [i in _types for i in road_types]
            ), f"`ROADS_kwargs['type']` must be a list of these: {_types}"
            road_geom = df.loc[
                df["type"].apply(lambda x: x in road_types)
            ].geometry.values
            self.ax.add_geometries(road_geom, crs=pc, **kwargs)

        if self.verbose == "debug":
            print("🐛 ROADS:", kwargs)
        return self

    def PLACES(
        self,
        country="United States",
        rank=2,
        scatter=True,
        labels=True,
        label_kw={},
        scatter_kw={},
    ):
        """
        Points and labels for major cities

        Parameters
        ----------
        country : str
            Country to filter
        rank : int
            City rank threshold. Large cities have small rank. Small
            cities have large rank.
        scatter : bool
            Add scatter points
        labels : bool
            Add city name labels
        """
        scatter_kw.setdefault("marker", ".")
        label_kw.setdefault("fontweight", "bold")
        label_kw.setdefault("alpha", 0.5)

        places = shapereader.natural_earth("10m", "cultural", "populated_places")
        df = geopandas.read_file(places)

        df = df[df.SOV0NAME == country]
        df = df[df.SCALERANK <= rank]

        xs = df.geometry.values.x
        ys = df.geometry.values.y
        names = df.NAME

        if scatter:
            self.ax.scatter(xs, ys, transform=pc, **scatter_kw)

        if labels:
            for x, y, name in zip(xs, ys, names):
                self.ax.text(x, y, name, clip_on=True, **label_kw)
        return self

    # Tiled images
    def STAMEN(self, style="terrain-background", zoom=3, alpha=1):
        """
        Add Stamen map tiles to background.

        .. note::
            When adding a tile product to a map, it might be better to add
            it to the map first, then set the map extent, then make a separate
            call to ``EasyMap()`` to add other features like roads and
            counties. The reason is because, if you add a tile map to

        Parameters
        ----------
        style : {'terrain-background', 'terrain', 'toner-background', 'toner', 'watercolor'}
            Type of image tile
        zoom : int
            Zoom level between 0 and 10.
        alpha : float
            Alpha value (transparency); a value between 0 and 1.
        """

        # Zoom can't be bigger than 11
        zoom = min(11, zoom)

        # Zoom can't be smaller than 0
        zoom = max(0, zoom)

        if self.verbose:
            print(
                "😎 Please use `ax.set_extent` before increasing Zoom level for faster plotting."
            )
        stamen_terrain = cimgt.Stamen(style)
        self.ax.add_image(stamen_terrain, zoom)

        if alpha < 1:
            # Need to manually put a white layer over the STAMEN terrain
            if self.dark:
                alpha_color = "k"
            else:
                alpha_color = "w"
            poly = self.ax.projection.domain
            self.ax.add_feature(
                feature.ShapelyFeature([poly], self.ax.projection),
                color=alpha_color,
                alpha=1 - alpha,
                zorder=1,
            )
        if self.verbose == "debug":
            print("🐛 STAMEN:", f"{style=}, {zoom=}, {alpha=}")

        return self

    def OSM(self, zoom=1, alpha=1):
        """
        Add Open Street Map tiles as background image.

        .. note::
            When adding a tile product to a map, it might be better to add
            it to the map first, then set the map extent, then make a separate
            call to ``EasyMap()`` to add other features like roads and
            counties. The reason is because, if you add a tile map to

        Parameters
        ----------
        zoom : int
            Zoom level between 0 and ?.
        alpha : float
            Alpha value (transparency); a value between 0 and 1.
        """

        image = cimgt.OSM()
        self.ax.add_image(image, zoom)
        if alpha < 1:
            # Need to manually put a white layer over the STAMEN terrain
            if self.dark:
                alpha_color = "k"
            else:
                alpha_color = "w"
            poly = self.ax.projection.domain
            self.ax.add_feature(
                feature.ShapelyFeature([poly], self.ax.projection),
                color=alpha_color,
                alpha=1 - alpha,
                zorder=1,
            )
        if self.verbose == "debug":
            print("🐛 OSM:", f"{zoom=}, {alpha=}")

        return self

    def STOCK(self, **kwargs):
        """Show stock image background (suitable for full-globe images)"""
        self.ax.stock_img()
        return self

    # Other
    def DOMAIN(
        self,
        x,
        y=None,
        *,
        text=None,
        method="cutout",
        facealpha=0.25,
        text_kwargs={},
        **kwargs,
    ):
        """
        Add a polygon of the domain boundary to a map.

        The border is drawn from the outside values of the latitude and
        longitude xarray coordinates or numpy array.
        Lat/lon values should be given as degrees.

        Parameters
        ----------
        x : xarray.Dataset or numpy.ndarray
            If xarray, then should contain 'latitude' and 'longitude' coordinate.
            If numpy, then 2D numpy array for longitude and `y` arg is required.
        y : numpy.ndarray
            Only required if x is a numpy array.
            A numpy array of latitude values.
        text : str
            If not None, puts the string in the bottom left.
        method : {'fill', 'cutout', 'border'}
            Plot the domain as a filled area Polygon, a Cutout from the
            map, or as a simple border.
        facealpha : float between 0 and 1
            Since there isn't a "facealpha" attribute for plotting,
            this will be it.
        polygon_only : bool
            - True: Only return the polygons and don't plot on axes.
        """
        _method = {"fill", "cutout", "border"}
        assert method in _method, f"Method must be one of {_method}."

        ####################################################################
        # Determine how to handle output...xarray or numpy
        if isinstance(x, (xr.core.dataset.Dataset, xr.core.dataarray.DataArray)):
            if self.verbose:
                print("process input as xarray")

            if "latitude" in x.coords:
                x = x.rename({"latitude": "lat", "longitude": "lon"})
            LON = x.lon.data
            LAT = x.lat.data

        elif isinstance(x, np.ndarray):
            assert y is not None, "Please supply a value for x and y"
            if self.verbose:
                print("process input as numpy array")
            LON = x
            LAT = y
        else:
            raise ValueError("Review your input")
        ####################################################################

        # Path of array outside border starting from the lower left corner
        # and going around the array counter-clockwise.
        outside = (
            list(zip(LON[0, :], LAT[0, :]))
            + list(zip(LON[:, -1], LAT[:, -1]))
            + list(zip(LON[-1, ::-1], LAT[-1, ::-1]))
            + list(zip(LON[::-1, 0], LAT[::-1, 0]))
        )
        outside = np.array(outside)

        ## Polygon in latlon coordinates
        ## -----------------------------
        x = outside[:, 0]
        y = outside[:, 1]
        domain_polygon_latlon = sgeom.Polygon(zip(x, y))

        ## Polygon in projection coordinates
        ## ----------------------------------
        transform = self.ax.projection.transform_points(pc, x, y)

        # Remove any points that run off the projection map (i.e., point's value is `inf`).
        transform = transform[~np.isinf(transform).any(axis=1)]

        # These are the x and y points we need to create the Polygon for
        x = transform[:, 0]
        y = transform[:, 1]

        domain_polygon = sgeom.Polygon(
            zip(x, y)
        )  # This is the boundary of the LAT/LON array supplied.
        global_polygon = (
            self.ax.projection.domain
        )  # This is the projection globe polygon
        cutout = global_polygon.difference(
            domain_polygon
        )  # This is the difference between the domain and glob polygon

        # Plot
        kwargs.setdefault("edgecolors", "k")
        kwargs.setdefault("linewidths", 1)
        if method == "fill":
            kwargs.setdefault("facecolor", (0, 0, 0, facealpha))
            artist = self.ax.add_feature(
                feature.ShapelyFeature([domain_polygon], self.ax.projection),
                **kwargs,
            )
        elif method == "cutout":
            kwargs.setdefault("facecolor", (0, 0, 0, facealpha))
            artist = self.ax.add_feature(
                feature.ShapelyFeature([cutout], self.ax.projection), **kwargs
            )
        elif method == "border":
            kwargs.setdefault("facecolor", "none")
            artist = self.ax.add_feature(
                feature.ShapelyFeature([domain_polygon.exterior], self.ax.projection),
                **kwargs,
            )

        if text:
            text_kwargs.setdefault("verticalalignment", "bottom")
            text_kwargs.setdefault("fontsize", 15)
            xx, yy = outside[0]
            self.ax.text(xx + 0.2, yy + 0.2, text, transform=pc, **text_kwargs)

        self.domain_polygon = domain_polygon
        self.domain_polygon_latlon = domain_polygon_latlon
        return self

    def INSET_GLOBE(self, **kwargs):
        """Add an axis showing a global inset map"""
        return inset_global_map(self.ax, **kwargs)

    def set_extent(self, *args, **kwargs):
        self.ax.set_extent(*args, **kwargs)
        return self

    def center_extent(self, *args, **kwargs):
        self.ax.center_extent(*args, **kwargs)
        return self

    def adjust_extent(self, *args, **kwargs):
        self.ax.adjust_extent(*args, **kwargs)
        return self

#===========================================================================
# MAIN PLOTTING FUNCTIONS
#===========================================================================

def plot_nwp(model, product, date_nwp, area, extent, fxx, var, scale, offset, level_hpa, level_min, level_max, level_int, plot_type, contour_color, apply_cmap, cmap,
             view_clabel, fontsize, linewidth, linestyle, land_ocean, land_color, ocean_color, alpha, title, figsize):

    # start the time counter
    start_counter = time.time()

    # convet date strings to ints
    year = int(date_nwp[0:4])
    month = int(date_nwp[5:7])
    day = int(date_nwp[8:10])
    hour = int(date_nwp[11:13])
    minute = int(date_nwp[14:16])

    # if using the forecast time as date, the plot will have the valid time as the date
    if fxx == 'date':
        multiple = 3
        number = hour
        # calculate the nearest fxx based in the hour
        fxx = multiple * round(number / multiple)
    else:
        # plot an specific forecast time (not the "valid" one from the date)
        fxx = fxx

    # create a Herbie object
    H = Herbie(f'{str(year).zfill(4)}-{str(month).zfill(2)}-{str(day).zfill(2)}', model=model, product=product, fxx=fxx)

    # check the desired variable and depending o variable, read the dataset accordingly
    # streamlines
    if var == 'winds':
        # read the dataset
        ds = H.xarray(f':(v|u):{level_hpa}')
        # calculate the windspeed
        ws = np.sqrt(ds.u**2 + ds.v**2)
        print(ws.shape)
    # wind speed
    elif var == 'windspeed':
        # read the dataset
        ds = H.xarray(f':(v|u):{level_hpa}')
        # calculate the windspeed
        ws = np.sqrt(ds.u**2 + ds.v**2)
    # thickness
    elif var == 'thickness':
        # read the dataset
        ds1 = H.xarray(f":gh:{level_hpa.split(',')[0]}")
        # read the dataset
        ds2 = H.xarray(f":gh:{level_hpa.split(',')[1]}")
    # vertical velocity
    elif var == 'VVEL':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 'w'
    # 2 m temperature
    elif var == '2t':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 't2m'
        if cmap == 'custom':
            # create a custom color scale:
            # The reference color palette may be found at the following page: http://wxmaps.org/pix/temp8
            # You may get the hex files st the following page: https://imagecolorpicker.com/
            colors = ["#d3d2d2", "#bcbcbc", "#969696", "#1464d2", "#1e6eeb", "#2882f0",
            "#3c96f5", "#50a5f5", "#78b9fa", "#96d2fa", "#b4f0fa", "#1eb41e", "#37d23c",
            "#50f050", "#78f573", "#96f58c", "#b4faaa", "#c8ffbe", "#ffe878", "#ffc03c",
            "#ffa000", "#ff6000", "#ff3200", "#e11400", "#c00000", "#a50000", "#785046",
            "#8c6359", "#b48b82", "#e1beb4"]
            cmap = matplotlib.colors.ListedColormap(colors)
            cmap.set_over('#fadad5')
            cmap.set_under('#e5e5e5')
    # low cloud cover
    elif var == 'LCDC':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 'lcc'
    # medium cloud cover
    elif var == 'MCDC':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 'mcc'
    # high cloud cover
    elif var == 'HCDC':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 'hcc'
    # total cloud cover
    elif var == 'TCDC':
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')
        var = 'tcc'
    # other variables
    else:
        # read the dataset
        ds = H.xarray(f':{var}:{level_hpa}')

    # declare ax as global (so it can be used be subsequent plots)
    global ax

    # check the desired area (previous or custom)
    if area == 'previous':
        print(f'\nUsing the same area of the previous plot...\n')
    elif area == 'custom':
        # choose the plot size (width x height, in inches)
        plt.figure(figsize=(figsize))
        # use the Cilindrical Equidistant projection in cartopy
        ax = plt.axes(projection=ccrs.PlateCarree())
        # define the image extent
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
        ax.set_extent([extent[0], extent[2], extent[1], extent[3]], crs=ccrs.PlateCarree())
        if land_ocean == True:
            # Add a land mask
            ax.add_feature(cfeature.LAND, facecolor=land_color)
            ax.add_feature(cfeature.OCEAN, facecolor=ocean_color)
        # add coastlines, borders and gridlines
        ax.coastlines(resolution='50m', color='black', linewidth=1.0)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1.0)
        # add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.5,
        xlocs=np.arange(-180, 181, 10), ylocs=np.arange(-90, 91, 10), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        gl.xpadding = -5
        gl.ypadding = -5
        gl.ylabel_style = {'color': 'white', 'size': 6, 'weight': 'bold'}
        gl.xlabel_style = {'color': 'white', 'size': 6, 'weight': 'bold'}

    # check the plot type. depending on the plot type, create the plot accordingly
    # unfilled contour
    if plot_type == 'c_unfilled':
        # thickness
        if var == 'thickness':
            if apply_cmap == False:
                # contours (unfilled) for ds
                img1 = ax.contour(ds1.longitude, ds1.latitude, ds2.gh/scale - ds1.gh*scale, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            elif apply_cmap == True:
                # contours (unfilled) for ds
                img1 = ax.contour(ds1.longitude, ds1.latitude, ds2.gh/scale - ds1.gh*scale, transform=pc, cmap=cmap, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img1, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
        # windspeed
        elif var == 'windspeed':
            if apply_cmap == False:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude, ws, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            elif apply_cmap == True:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude, ws, transform=pc, cmap=cmap, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img1, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color, alpha=alpha)
        # low, medium or high cloud cover
        elif var == 'lcc' or var == 'mcc' or var == 'hcc':
            if apply_cmap == False:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude, (getattr(ds, var.lower())*scale+offset)[0], transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            elif apply_cmap == True:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude,  (getattr(ds, var.lower())*scale+offset)[0], transform=pc, cmap=cmap, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img1, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
        # other variables
        else:
            if apply_cmap == False:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude, getattr(ds, var.lower())*scale+offset, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            elif apply_cmap == True:
                # contours (unfilled) for ds
                img1 = ax.contour(ds.longitude, ds.latitude, getattr(ds, var.lower())*scale+offset, transform=pc, cmap=cmap, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img1, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)

    # filled contour
    elif plot_type == 'c_filled':
        # thickness
        if var == 'thickness':
            # contours (filled) for ds
            img1 = ax.contourf(ds.longitude, ds.latitude, ds2.gh*scale - ds1.gh*scale, transform=pc, cmap=cmap, levels=np.arange(level_min, level_max, level_int), extend='both', alpha=alpha)
            # contours (unfilled) for ds
            img2 = ax.contour(ds.longitude, ds.latitude, ds2.gh*scale - ds1.gh*scale, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img2, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
        # windspeed
        elif var == 'windspeed':
            # contours (filled) for ds
            img1 = ax.contourf(ds.longitude, ds.latitude, ws, transform=pc, cmap=cmap, levels=np.arange(level_min, level_max, level_int), extend='both', alpha=alpha)
            # contours (unfilled) for ds
            img2 = ax.contour(ds.longitude, ds.latitude, ws, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img2, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
        # low, medium or high cloud cover
        elif var == 'lcc' or var == 'mcc' or var == 'hcc':
            # contours (filled) for ds
            img1 = ax.contourf(ds.longitude, ds.latitude,  (getattr(ds, var.lower())*scale+offset)[0], transform=pc, cmap=cmap, levels=np.arange(level_min, level_max, level_int), extend='both', alpha=alpha)
            if view_clabel == True:
                # labels for ds
                ax.clabel(img1, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
        # other variables
        else:
            # contours (filled) for ds
            img1 = ax.contourf(ds.longitude, ds.latitude, getattr(ds, var.lower())*scale+offset, transform=pc, cmap=cmap, levels=np.arange(level_min, level_max, level_int), extend='both', alpha=alpha)
            # contours (unfilled) for ds
            img2 = ax.contour(ds.longitude, ds.latitude, getattr(ds, var.lower())*scale+offset, transform=pc, colors=contour_color, linewidths=linewidth, linestyles=linestyle, levels=np.arange(level_min, level_max, level_int), alpha=alpha)
            
            if view_clabel == True:
                # labels for ds
                ax.clabel(img2, inline=1, inline_spacing=0, fontsize=fontsize, fmt = '%1.0f', colors=contour_color)
    
    # pcolormesh
    elif plot_type == 'pmesh':
      # pcolormesh for ds
      p = ax.pcolormesh(ds.longitude, ds.latitude, getattr(ds, var.lower())*scale+offset, transform=pc, cmap=cmap, vmin=level_min, vmax = level_max, alpha=alpha)
  
    # streamlines
    elif plot_type == 'streamlines':
      #if apply_cmap == True:
        # streamlines for ds
      #  img1 = ax.streamplot(ds.longitude, ds.latitude, ds.u*scale, ds.v*scale, density=[4, 4], color=ws, cmap=cmap, linewidth=linewidth, transform=pc)
      #else:
        # streamlines for ds
        img1 = ax.streamplot(ds.longitude, ds.latitude, ds.u*scale, ds.v*scale, density=[4, 4], linewidth=linewidth, transform=pc, color=contour_color)
  
    # wind barbs
    elif plot_type == 'barbs':
      if apply_cmap == True:
         # barbs for ds
        flip_flag = np.zeros((ds.u.shape[0],ds.u.shape[1]))
        flip_flag[ds.latitude < 0] = 1
        img1 = ax.barbs(ds.longitude[::8], ds.latitude[::8], ds.u[::8,::8], ds.v[::8,::8], ws[::8,::8], cmap=cmap, transform=pc, length = 4.0, sizes = dict(emptybarb=0.0, spacing=0.2, height=0.5), linewidth=linewidth, linestyles=linestyle, pivot='middle', flip_barb = flip_flag[::8,::8], alpha=alpha)
      else:
        # barbs for ds
        flip_flag = np.zeros((ds.u.shape[0],ds.u.shape[1]))
        flip_flag[ds.latitude < 0] = 1
        img1 = ax.barbs(ds.longitude[::8], ds.latitude[::8], ds.u[::8,::8], ds.v[::8,::8], transform=pc, length = 4.0, sizes = dict(emptybarb=0.0, spacing=0.2, height=0.5), linewidth=linewidth, linestyles=linestyle, pivot='middle', barbcolor=contour_color, flip_barb = flip_flag[::8,::8], alpha=alpha)
  
    # wind vectors
    elif plot_type == 'vectors':
      if apply_cmap == True:
        # vectors for ds
        img1 = ax.quiver(ds.longitude[::8], ds.latitude[::8], ds.u[::8,::8], ds.v[::8,::8],  ws[::8,::8], cmap=cmap)
      else:
        # vectors for ds
        img1 = ax.quiver(ds.longitude[::8], ds.latitude[::8], ds.u[::8,::8], ds.v[::8,::8], color=contour_color)
  
    # check if we want to concatenate the variable name (multiple plot) or not
    if title == 'new':
      global variable_title
      variable_title = ''
  
    # add title and other decorations
    if var == 'thickness':
      variable_title = variable_title + ' | ' + f'Thickness {ds1.isobaricInhPa.item()} - {ds2.isobaricInhPa.item()} {ds1.isobaricInhPa.units}'
      plt.title(f"{ds1.model.upper()}: {H.product_description}\nValid: {ds1.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # mean sea level
    elif var == 'msl':
      # use definition to plot H/L symbols
      plot_maxmin_points(ds.longitude[::2], ds.latitude[::2], ds.msl[::2,::2]/100, 'max', 50, symbol='H', color='b',  transform=pc)
      plot_maxmin_points(ds.longitude[::2], ds.latitude[::2], ds.msl[::2,::2]/100, 'min', 25, symbol='L', color='r', transform=pc)
      variable_title = variable_title + ' | ' + f'{getattr(ds, var.lower()).GRIB_name}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # winds or streamlines
    elif var == 'winds' and plot_type == 'streamlines':
      variable_title = variable_title + ' | ' +  f'Streamlines {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # windspeed
    elif var == 'windspeed':
      variable_title = variable_title + ' | ' +  f'Windspeed (m s⁻¹) {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # barbs
    elif var == 'winds' and plot_type == 'barbs':
      variable_title = variable_title + ' | ' +  f'Wind Barbs {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
     # barbs
    elif var == 'winds' and plot_type == 'vectors':
      variable_title = variable_title + ' | ' +  f'Wind Vectors {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # total precipitation
    elif var == 'tp':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} (mm)'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # 2 m temperature
    elif var == 't2m':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} ({getattr(ds, var.lower()).units})'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # temperature
    elif var == 't':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} (°C) {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # convective precipitation
    elif var == 'ACPCP':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} ({getattr(ds, var.lower()).units})'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # low, medium or high cloud cover
    elif var == 'lcc' or var == 'mcc' or var == 'hcc':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} ({getattr(ds, var.lower()).units})'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # precipitable water
    elif var == 'PWAT':
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} ({getattr(ds, var.lower()).units})'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
    # other variables
    else:
      variable_title = variable_title + ' | ' +  f'{getattr(ds, var.lower()).GRIB_name} ({getattr(ds, var.lower()).units}) {ds.isobaricInhPa.item()} {ds.isobaricInhPa.units}'
      plt.title(f"{ds.model.upper()}: {H.product_description}\nValid: {ds.valid_time.dt.strftime('%H:%M UTC %d %b %Y').item()}\n{variable_title}",loc='right')
  
    # save the image
    img_file = f'output/{model}_{var}_{date_nwp.replace(" ", "_")}.png'
    plt.savefig(f'{img_file}', bbox_inches='tight', pad_inches=0, dpi=100)
  
    # show the image
    #plt.show()
  
    # print to total time required to process the file
    print(f'Total processing time (nwp): {round((time.time() - start_counter),2)} seconds.')
    print(f'--------------------------------------------------------------------------------------------------------------------------------------')
  
    # return the image file name
    return img_file

#===========================================================================
# MAXIMUM/MINIMUM POINT PLOTTING
#===========================================================================

def plot_maxmin_points(lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
    """
    This function will find and plot relative maximum and minimum for a 2D grid. The function
    can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
    values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
    max/min value. The symbol text can be set to a string value and optionally the color of the
    symbol and any plotted value can be set with the parameter color
    lon = plotting longitude values (2D)
    lat = plotting latitude values (2D)
    data = 2D data that you wish to plot the max/min symbol placement
    extrema = Either a value of max for Maximum Values or min for Minimum Values
    nsize = Size of the grid box to filter the max and min values to plot a reasonable number
    symbol = String to be placed at location of max/min value
    color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
    plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
    The max/min symbol will be plotted on the current axes within the bounding frame
    (e.g., clip_on=True)
    """
    # import scipy filters
    from scipy.ndimage.filters import maximum_filter, minimum_filter

    # algorithm to find minimuns and maximuns
    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    # plot the high and low labels
    for i in range(len(mxy)):
         txt1 = ax.annotate(symbol, xy=(lon[mxx[i]], lat[mxy[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), color=color, size=16,
                fontweight='bold', horizontalalignment='center', verticalalignment='center', transform=pc, annotation_clip=True, clip_on=True)

         txt2 = ax.annotate('\n' + str(int(data[mxy[i], mxx[i]])), xy=(lon[mxx[i]], lat[mxy[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                color=color, size=8, fontweight='bold', horizontalalignment='center', verticalalignment='top', transform=pc, annotation_clip=True, clip_on=True)

#===========================================================================
# SCRIPT EXECUTION - CUSTOMIZE PARAMETERS HERE
#===========================================================================

plot_nwp(model='ecmwf',
         product='oper',
         date_nwp='2024-02-13 12:00',
         area='custom',
         extent=[-20.0, 30.00, 25.00, 60.00],
         fxx='date',
         var='msl', scale=0.01, offset=0,
         level_hpa='',
         level_min=990, level_max=1050, level_int=4,
         plot_type='c_unfilled', contour_color='white',
         apply_cmap=False, cmap='jet',
         view_clabel=True, fontsize=6,
         linewidth=0.5, linestyle='solid',
         land_ocean=True, land_color = 'lightgray', ocean_color = 'white',
         alpha=1.0, title='new',
         figsize=[15,15])