"""
Created 22 Jan 2021
Last updated 12 Feb 2021
GIS functions for kartturs GeoImagine Framework

Author
______
Thomas Gumbricht
"""

from .version import __version__, VERSION, metadataD

from .kt_gis import MjProj, GetVectorProjection, GetRasterMetaData, Geometry, ESRIOpenGetLayer, RasterOpenGetFirstLayer

__all__ = ['MjProj','GetVectorProjection','GetRasterMetaData','Geometry','ESRIOpenGetLayer']