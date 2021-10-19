'''
Created on 21 Jan 2021
Last updated 12 Feb 2021

@author: thomasgumbricht
'''

# Standard library imports

from __future__ import division

import os

import sys

import array as arr

from itertools import chain

# Third party imports

import numpy as np

from osgeo import gdal, ogr, osr

from gdalconst import *

#import fiona
#import json

# Package application imports

#import mj_datetime_v70 as mj_dt

class MjProj:
    ''' Class for handling projections
    '''
    
    def __init__(self):
        """The constructor is just an empty container.""" 
        
        pass
         
    def SetProj(self,spatialRef):
        '''
        '''
        
        self.proj_cs = spatialRef
        
    def SetGeoTrans(self,gt):
        '''
        '''
        self.gt = gt
        
    def SetFromProj4(self,proj4):
        '''
        '''
        
        self.proj4 = proj4
        
        srs = osr.SpatialReference()
          
        srs.ImportFromProj4(self.proj4)
        
        self.proj_cs = osr.SpatialReference()
        
        self.proj_cs.ImportFromWkt(srs.ExportToWkt())
        
    def SetProj4(self): 
        '''
        '''
        
        srs = osr.SpatialReference()
        
        self.WKT = self.proj_cs.ExportToWkt()
        
        srs.ImportFromWkt(self.WKT)
        
        self.proj4 = srs.ExportToProj4()
        
    def SetFromWKT(self,WKT):
        '''
        '''
        
        self.WKT = WKT
        
        self.proj_cs = osr.SpatialReference()
        
        self.proj_cs.ImportFromWkt(WKT)
        
    def SetFromEPSG(self,epsg):
        '''
        '''
        
        self.epsg = epsg 
             
        self.proj_cs = osr.SpatialReference() 
        
        if epsg == 6842:
            
            WKT = 'PROJCS["Sinusoidal_Sanson_Flamsteed",GEOGCS["GCS_Unknown",DATUM["D_unknown",SPHEROID["Unknown",6371007.181,"inf"]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Sinusoidal"],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
           
            self.SetFromWKT(WKT)
        
        else:
            
            self.proj_cs.ImportFromEPSG(epsg)
            
    def ReadSpatialRef(self):
        '''
        '''
        self.SetProj4()

        self.SetFromProj4(self.proj4)
        if self.proj_cs.GetAttrValue("AUTHORITY", 1) == None:
            self.epsg = False
            return
        self.epsg = int( self.proj_cs.GetAttrValue("AUTHORITY", 1) )
        self.spatialunit = self.proj_cs.GetAttrValue("UNIT", 0)
        self.datum = self.proj_cs.GetAttrValue("DATUM", 0)
        self.spheroid = self.proj_cs.GetAttrValue("SPHEROID", 0)
        self.projection = self.proj_cs.GetAttrValue("PROJECTION", 0)
        self.setEpsg = False
        if self.epsg == None and self.datum in ['WGS_1984','WGS_84'] and self.spheroid in ['WGS_1984','WGS_84'] and self.spatialunit.lower() == 'degree':
            self.epsg = 4326
            self.setEpsg = True
        elif self.epsg == 4326:
            pass
        elif self.proj4.replace(' ','') == '+proj=sinu+lon_0=0+x_0=0+y_0=0+a=6371007.181+b=6371007.181+units=m+no_defs':
            self.epsg = 6842
            self.setEpsg = True
        elif self.projection == None:
            self.epsg = False
 
class ogrFeature:
    def __init__(self,layer):
        """The constructoris just an empty container.""" 
        self.layer = layer.layer
         
    def SetFeature(self,feature):
        self.feature = feature 
        
    def CreateOgrFeature(self, geom, fieldDefL, dictL = False, key = False):
        '''
        '''
        featDefn = self.layer.GetLayerDefn()
        tarfeat = ogr.Feature(featDefn)
        tarfeat.SetGeometry(geom)
        
        for f in fieldDefL:
            if f.transfer in ['constant','fixed','field']:
                tarfeat.SetField(f.name,f.source)       
            elif f.transfer == 'copy':
                exit( 'alternative copy not identified in kt_gis.ogrFeature.CreateOgrFeature')
            elif f.transfer == 'db':
                tarfeat.SetField(f.name,f.source)
            elif f.transfer == 'dict':
                tarfeat.SetField(f.name,f.source[key])
            else:
                exitstr = 'Unrecognized transfer command: %s' %(f.transfer)
                sys.exit(exitstr)
        self.layer.CreateFeature(tarfeat)
        
    def CopyOgrFeature(self,srcfeat,fieldDefL):
        featDefn = self.layer.GetLayerDefn()
        tarfeat = ogr.Feature(featDefn)
        tarfeat.SetGeometry(srcfeat.GetGeometryRef())
        for f in fieldDefL:        
            tarfeat.SetField(f.name, f.source)           
        self.layer.CreateFeature(tarfeat)      
        
    def AddSourceValtoFieldDef(self,fieldDefL):
        for f in fieldDefL:
            if f.transfer == 'fixed':
                pass
            elif f.transfer == 'field':

                f.source = self.feature.GetField(f.field)
            else:
                f.source = self.feature.GetField(f.name)
            
class Geometry:
    '''
    '''
    
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def GeomFromFeature(self,feature): 
        self.ogrGeom = feature.GetGeometryRef() 
        self.OgrGeomToShapely()  
           
    def SetShapelyGeom(self,geom):   
        self.shapelyGeom = geom
        self.ShapelyToOgrGeom()
        
    def LoadWKTGeom(self,geom):
        import shapely.wkt
        self.shapelyGeom = shapely.wkt.loads(geom)
        self.ShapelyToOgrGeom()
        
    def BoundsToPoly(self):
        from shapely.geometry import Polygon
        minx, miny, maxx, maxy = self.shapelyGeom.bounds
        ptLT = ( (minx,maxy), (maxx,maxy), (maxx,miny), (minx,miny), (minx,maxy) )
        return Polygon(ptLT)
    
    def PointsToMultiPointGeom(self,ptL):
        from shapely.geometry import Point,MultiPoint
        ptgeomL = []
        for pt in ptL:
            ptgeomL.append(Point(pt))
        self.shapelyGeom = MultiPoint(ptgeomL)
        self.ShapelyToOgrGeom()
        
    def PointsToPolygonGeom(self,ptL):
        from shapely.geometry import Polygon
        #ptgeomL = []
        #for pt in ptL:
        #    ptgeomL.append(Point(pt))
        self.shapelyGeom = Polygon(ptL)
        #self.shapelyGeom = Polygon(ptgeomL)
        self.ShapelyToOgrGeom()
        
    def CascadeUnion(self,shapelyGeomL):
        from shapely.ops import cascaded_union
        self.SetShapelyGeom(cascaded_union(shapelyGeomL))
        #return cascaded_union(shapelyGeomL)
        
    def MultPolyToSinglePoly(self,multipoly):
        from shapely.ops import cascaded_union
        from shapely.geometry import Polygon
        '''
        self.shapelyGeom = cascaded_union([
            Polygon(component.exterior) for component in multipoly
        ])
        '''
        polyL = [item for item in multipoly]
        self.shapelyGeom = polyL[0]
        for item in polyL:
            self.shapelyGeom.union(item)
        self.ShapelyToOgrGeom()
        
    def OgrGeomToShapely(self):
        from shapely.wkb import loads
        self.shapelyGeom = loads(self.ogrGeom.ExportToWkb())
        
    def ShapelyToOgrGeom(self):
        self.ogrGeom = ogr.CreateGeometryFromWkb(self.shapelyGeom.wkb)
        
    def ShapelyBuffer(self,buff):
        if not self.shapelyGeom.is_valid:
            print ( '    mending invalid geometry by buffer = 0' )
            self.shapelyGeom = self.shapelyGeom.buffer(0)
            if not self.shapelyGeom.is_valid:
                exitstr = 'Can not fix invalid geometry in ShapelyBuffer'
                sys.exit(exitstr)
                
    def ShapelyPolysToMultiPoly(self,multiPolyL):
        from shapely.geometry.multipolygon import MultiPolygon
        self.shapelyGeom = MultiPolygon(multiPolyL) 
        
    def ShapelyIntersection(self, otherGeom): 
        return self.shapelyGeom.intersection(otherGeom.shapelyGeom)
    
        
    def MultiPolyGeomFromGeomL(self,geomL):
        singlegeomL = []
        for testgeom in geomL:
            if testgeom.geom_type == 'MultiPolygon':
                for polygeom in [polygeom for polygeom in testgeom]:
                    if len(list(polygeom.exterior.coords)) > 2:
                        singlegeomL.append(polygeom)
            else:
                if len(list(testgeom.exterior.coords)) > 2:
                    singlegeomL.append(testgeom)
        self.ShapelyPolysToMultiPoly(self,singlegeomL)
        
    def SplitPolygonsToSingleGeomL(self):
        singlegeomL = []
        if self.shapelyGeom.geom_type == 'MultiPolygon':
            for polygeom in [polygeom for polygeom in self.shapelyGeom]:
                if len(list(polygeom.exterior.coords)) > 2:
                    singlegeomL.append(polygeom)
        else:
            if len(list(self.shapelyGeom.exterior.coords)) > 2:
                    singlegeomL.append(self.shapelyGeom)
        return singlegeomL
             
    def ShapelyRemoveSlivers(self):
        eps = 0.000001
        from shapely.geometry import JOIN_STYLE
        self.ShapelyBuffer(0)
        #remove slivers
        self.shapelyGeom = self.shapelyGeom.buffer(eps, 1, join_style=JOIN_STYLE.mitre).buffer(-eps, 1, join_style=JOIN_STYLE.mitre)
                
    def ShapelyContainCut(self,cut):
        if cut.shapelyGeom.contains(self.shapelyGeom):
            pass
        #geom remains intact
        elif cut.shapelyGeom.intersects(self.shapelyGeom):
            self.shapelyGeom = cut.shapelyGeom.intersection(self.shapelyGeom)
        else:
            self.shapelyGeom = False
                
    def ShapelyOutsideCut(self,cut):
        if cut.shapelyGeom.contains(self.shapelyGeom):
            self.shapelyGeom = False
        elif cut.shapelyGeom.intersects(self.shapelyGeom):
            self.shapelyGeom = self.shapelyGeom.difference(cut.shapelyGeom)
        else:
            pass
        #geom remains intact
        
    def ShapelyWithin(self,cut):
        if self.shapelyGeom.within(cut.shapelyGeom):
            return True
        else:
            return False
        
    def GeoTransform(self,srcProj,tarProj):
        self.geoTransform = osr.CoordinateTransformation(srcProj.proj_cs,tarProj.proj_cs)
        res = self.ogrGeom.Transform(self.geoTransform)
        if res != 0:
            sys.exit('Geotransformation failed')
        #set shapley to reproejcted
        self.OgrGeomToShapely()
            
    def GeoTransformCoordsFIX(self,srcProj,tarProj,ptL): 
        #Special version that allows geotransformation outised lat lon boundaries (i.e. for SIN outised erait sphere
        transform = osr.CoordinateTransformation(srcProj.proj_cs,tarProj.proj_cs)
        xyL = []
        for pt in ptL:
            xyL.append(transform.TransformPoint(pt[0],pt[1]))
        return xyL
    
class ShapelyPointGeom(Geometry):
    def __init__(self,pt):
        from shapely.geometry import Point
        Geometry.__init__(self)
        self.shapelyGeom = Point(pt)
        
class ShapelyMultiPointGeom(Geometry):
    def __init__(self,ptL):
        from shapely.geometry import MultiPoint
        Geometry.__init__(self)
        self.shapelyGeom = MultiPoint(ptL)

class ShapelyLineGeom(Geometry):
    def __init__(self,ptLT):
        from shapely.geometry import LineString
        Geometry.__init__(self)
        self.shapelyGeom = LineString(ptLT)

class ShapelyPolyGeom(Geometry):
    def __init__(self,ptLT):
        from shapely.geometry import Polygon
        Geometry.__init__(self)
        self.shapelyGeom = Polygon(ptLT)       
        
class FieldDef:
    def __init__(self, name,fieldDefD):
        '''{'type','transfer','source','width','precision','keyfield/field'}'''
        self.name = name
        if 'type' in fieldDefD:
            self.type = fieldDefD['type']
        if 'transfer' in fieldDefD:
            self.transfer = fieldDefD['transfer']
        if 'source' in fieldDefD:
            self.source = fieldDefD['source']
        if 'width' in fieldDefD:
            self.width = fieldDefD['width']
        else:
            self.width = 8
        if 'precision' in fieldDefD:
            self.precision = fieldDefD['precision']
        else:
            self.precision = 0
        if 'keyfield' in fieldDefD:
            self.keyfield = fieldDefD['keyfield']
        if 'field' in fieldDefD:
            self.keyfield = fieldDefD['field']
        else:
            self.keyfield = False
        
class VectorLayer:
    def __init__(self):
        """The constructor is just an empty container.""" 
        
    #VectorDataSource.__init(self)   
    def SetDS(self,DS): 
        self.datasource =  DS.datasource
             
    def GetLayer(self):
        geomTypeD = {1:'point',2:'line',3:'polygon'} 
        self.layer = self.datasource.GetLayer()
        self.geomtype = geomTypeD[self.layer.GetGeomType()]
        self.layerid = self.layer.GetName()
        self.GetLayerDef()
        
    def GetLayerDef(self):
        layerDefn = self.layer.GetLayerDefn()
        fieldDefD = {}  
        self.fieldDefL = []
        for i in range( layerDefn.GetFieldCount() ):            
            n = layerDefn.GetFieldDefn(i).GetName()
            t = layerDefn.GetFieldDefn(i).GetType()
            #fieldDefD['typename'] = layerDefn.GetFieldDefn(i).GetFieldTypeName()
            t = layerDefn.GetFieldDefn(i).GetFieldTypeName(t)
            w = layerDefn.GetFieldDefn(i).GetWidth()
            p = layerDefn.GetFieldDefn(i).GetPrecision()
            fieldDefD[n] = {'type':t.lower(), 'width':w, 'precision': p, 'transfer': 'constant'}      
            #fieldDefD['source'] = self.feature.GetField(fieldDefD['name'])
        for key in fieldDefD:
            self.fieldDefL.append(FieldDef(key,fieldDefD[key]))
                   
    def CreateOgrLayer(self,geomtype,tarLayerId):
        self.geomtype = geomtype
        if self.geomtype.lower() == 'point':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPoint) 
        elif self.geomtype.lower() == 'line':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbLineString) 
        elif self.geomtype.lower() == 'polygon':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPolygon) 
        else:
            exitstr = 'Can not recoginze geometry type %s for creating layer' %(self.geomtype) 
            sys.exit(exitstr)
        if self.layer is None:
            sys.exit('Could not create layer') 
            
    def AddFieldDef(self,fieldDef):
        #Get all field
        self.layerDefn = self.layer.GetLayerDefn()
        if fieldDef.name not in self.schemaD:
            #Check if field exists:

            if fieldDef.type == 'string':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTString)
                new_field.SetWidth(fieldDef.width)
            elif fieldDef.type == 'integer':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTInteger)
            elif fieldDef.type == 'real':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTReal)
            elif fieldDef.type == 'date':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTDate)
            elif fieldDef.type == 'time':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTTime)
            else:
                printstr = ('field type not recognized %s' %(fieldDef.type) )
                sys.exit(printstr)
            self.layer.CreateField(new_field)
        else:
            for i in range( self.layerDefn.GetFieldCount() ):
                
                if self.layerDefn.GetFieldDefn(i).GetName() == fieldDef.name:
                    fieldTypeCode = self.layerDefn.GetFieldDefn(i).GetType()
                    fieldType = self.layerDefn.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
                    fieldWidth = self.layerDefn.GetFieldDefn(i).GetWidth()
                    fieldPrecision = self.layerDefn.GetFieldDefn(i).GetPrecision()

                    if fieldType.lower() != fieldDef.type.lower():
                        sys.exit('Existing field of different type already exists')
         
    def AddField(self,f):
        pass
              
    def FillAllFields(self,fieldDefL,keyfield):
        for f in fieldDefL:
            if f.transfer == 'constant':
                for feature in self.layer:  
                    feature.SetField(f.name, f.source)
                    self.layer.SetFeature(feature)
            elif f.transfer == 'copy':
                #all features have the same value
                for feature in self.layer:  
                    copy = feature.GetField(f.source)
                    feature.SetField(f.name, copy)
                    self.layer.SetFeature(feature)
            elif f.transfer == 'dict':
                #all features have the same value
                for feature in self.layer:  
                    featkey = feature.GetField(keyfield)
                    fielddata = f.source[featkey]
                    feature.SetField(f.name, fielddata)
                    self.layer.SetFeature(feature)
            else:
                exitstr = 'missing transfer type %s' %(f.transfer)
                sys.exit(exitstr)
       
    def AddPtDataFromDict(self,xkey,ykey,dictL):
        '''
        '''

        featDefn = self.layer.GetLayerDefn()
        
        for D in dictL:
            
            tarfeat = ogr.Feature(featDefn)
            
            #set the geometry
            
            geom = ShapelyPointGeom( (D[xkey], D[ykey]) )
            
            geom.ShapelyToOgrGeom()
            
            tarfeat.SetGeometry(geom.ogrGeom)   
               
            for key in D:  
                
                tarfeat.SetField(key, D[key]) 
                        
            self.layer.CreateFeature(tarfeat) 
       
    def _GetSpatialRef(self): 
        self.spatialRef = self.layer.GetSpatialRef() 
        
    def SetSpatialRef(self, tarProj):
        self.spatialRef = tarProj.proj_cs
        
    def SetCSSpatialRef(self, proj_cs):
        self.spatialRef = proj_cs
        
    def GetLayerFields(self):  
        self.schemaD = {}
        ldefn = self.layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            self.schemaD[fdefn.name]= [fdefn.GetType(),fdefn.GetWidth(),fdefn.GetPrecision()]
        
    def CompareFields(self,srcFieldId,srcFieldName,fieldDefL):
        res = True
        for feature in self.layer:  
            item = feature.GetField(srcFieldId)
            for f in fieldDefL:
                if f.transfer == 'dict':
                    if item in f.source:
                        pass
                    else:
                        print ('missing', item,feature.GetField(srcFieldName))
                        res = False
        return res  
               
class VectorDataSource:  
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def OpenESRIshapeRead(self,FPN): 
        self.ESRIshapeFPN = FPN
        if os.path.exists(self.ESRIshapeFPN):
            driver = ogr.GetDriverByName("ESRI Shapefile")
            self.datasource = driver.Open(self.ESRIshapeFPN, 0)
        else:
            exitstr = 'ESRI datasource %s does not exist' %(self.ESRIshapeFPN)
            sys.exit(exitstr)
            
    def OpenESRIshapeEdit(self,FPN): 
        self.ESRIshapeFPN = FPN
        if os.path.exists(self.ESRIshapeFPN):
            driver = ogr.GetDriverByName("ESRI Shapefile")
            self.datasource = driver.Open(self.ESRIshapeFPN, 1)
        else:
            exitstr = 'ESRI datasource %s does not exist' %(self.ESRIshapeFPN)
            sys.exit(exitstr)
            
    def CreateESRIdatasource(self,FPN):
        self.FPN = FPN
        driver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(self.FPN):
            driver.DeleteDataSource(self.FPN)
        self.datasource = driver.CreateDataSource(self.FPN)
        if self.datasource is None:
            print ('Could not create file')
            sys.exit(1)
                   
    def CloseDS(self):
        self.datasource = None 
        

        
class RasterDataSource: 
    def __init__(self):
        """The constructoris just an empty container.""" 
        pass
        
    def OpenGDALRead(self,FPN): 
        self.rasterFPN = FPN
        if os.path.exists(self.rasterFPN):
            self.datasource = gdal.Open(FPN, GA_ReadOnly)
            if self.datasource is None:
                exitstr = 'Exiting - Failed to open raster file %s' %(self.rasterFPN)
                sys.exit(exitstr)
        else:
            exitstr = 'Raster datasource %s does not exist' %(self.rasterFPN)
            sys.exit(exitstr)
                                                      
    def OpenGDALEdit(self,FPN): 
        self.rasterFPN = FPN
        if os.path.exists(self.rasterFPN):
            self.datasource = gdal.Open(FPN)
            if self.datasource is None:
                exitstr = 'Exiting - Failed to open raster file %s' %(self.rasterFPN)
                sys.exit(exitstr)
        else:
            exitstr = 'Raster datasource %s does not exist' %(self.rasterFPN)
            sys.exit(exitstr)
       
    
        
    def CreateGDALDS(self,FPN,layer,of = 'GTiff'):
        #self.datasource = gdal.GetDriverByName('GTiff').Create(FPN, x_res, y_res, 1, gdal.GDT_Byte)
        
        ''' NOTE that layerl.lins, layer.cols changed to layer.layer.cols and layer.layer.lins 20210329
        '''
        driver = gdal.GetDriverByName( of )
        metadata = driver.GetMetadata()
        #if metadata.has_key(gdal.DCAP_CREATE) \
        if gdal.DCAP_CREATE in metadata \
            and metadata[gdal.DCAP_CREATE] == 'YES':
            if layer.comp.celltype.lower() in ['byte','uint8']:

                self.datasource = driver.Create( FPN, layer.layer.cols, layer.layer.lins, 1, gdal.GDT_Byte )   
                layer.comp.cellnull = int(layer.comp.cellnull) 
            elif layer.comp.celltype.lower() == 'int16':  
                self.datasource = driver.Create( FPN, layer.layer.cols, layer.layer.lins, 1, gdal.GDT_Int16 ) 
                layer.comp.cellnull = int(layer.comp.cellnull)
            elif layer.comp.celltype.lower() == 'uint16':
                self.datasource = driver.Create( FPN, layer.layer.cols, layer.layer.lins, 1, gdal.GDT_UInt16 )
                layer.comp.cellnull = int(layer.comp.cellnull)
            elif layer.comp.celltype.lower() == 'float32':
                self.datasource = driver.Create( FPN, layer.layer.cols, layer.layer.lins, 1, gdal.GDT_Float32 )
                layer.comp.cellnull = float(layer.comp.cellnull)
            else:
                exitstr ='numpy type not defined',layer.comp.celltype
                sys.exit(exitstr) 
            self.datasource.SetGeoTransform( layer.layer.geotrans )
            self.datasource.SetProjection( layer.layer.projection  )
            self.datasource.GetRasterBand(1).SetNoDataValue(layer.comp.cellnull)
            if hasattr(layer.comp,'palette'):
                #PcR,AT,maxAT = self.FixGDALPalette(self.comp.palette)
                #ct = gdal.ColorTable() 
                #for c in range(1,len(PcR)):
                #    ct.CreateColorRamp(PcR[c-1][0],PcR[c-1][1],PcR[c][0],PcR[c][1])
                
                self.datasource.GetRasterBand(1).SetColorTable(layer.comp.palette.colortable)
                
    def _WriteFullArray(self,layer):
        '''
        '''
        # Write the data array
        self.datasource.GetRasterBand(1).WriteArray( layer.layer.NPBAND )
        
        # Recalculate the stats
        self._SetStats()
                        
    def _SetStats(self):
        ''' Calculate statistics
        '''
        # Get the stats
        stats = self.datasource.GetRasterBand(1).ComputeStatistics(0)
           
        # Write the stats     
        self.datasource.GetRasterBand(1).SetStatistics(stats[0],stats[1],stats[2],stats[3])
                    
    def WriteBlock(self,col,row,arr):
        self.datasource.GetRasterBand(1).WriteArray(arr, col, row)
        #self.datasource.GetRasterBand(1).WriteArray
 
    def arrTo2Dnp(self,layer):
        if layer.comp.celltype.lower() in ['byte','uint8']:
            band2DArray = np.asarray(layer.BAND, dtype=np.int8)
        elif layer.comp.celltype.lower() == 'int16':  
            band2DArray = np.asarray(layer.BAND, dtype=np.int16)    
        elif layer.comp.celltype.lower() == 'uint16':
            band2DArray = np.asarray(layer.BAND, dtype=np.uint16)
        else:
            ( 'print numpy type not defined',layer.comp.celltype )
            sys.exit()   
        #reshape the 1D array to a 2D image
        band2DArray.shape = (-1, layer.cols)
        layer.BAND = band2DArray
    
    def SetGeoTransform(self,gt): 
        self.datasource.SetGeoTransform(gt) 
    
    def SetProjection(self,proj): 
        self.datasource.SetProjection(proj)  
            
    def CloseDS(self):
        #close the datasource
        self.datasource.FlushCache()
        self.datasource = None 
        
class RasterLayer:
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def SetDS(self,DS): 
        self.datasource =  DS.datasource
        
    def GetLayer(self,bandnr):
        self.layer = self.datasource.GetRasterBand(bandnr)
        
    def GetSpatialRef(self): 
        self.spatialRef = self.proj_cs = self.datasource.GetProjection()
        self.gt = self.datasource.GetGeoTransform()
        
        #self.spatialRef = self.proj_cs = self.layer.GetSpatialRef() 
    def SetSpatialRef(self, tarProj):
        self.spatialRef = tarProj.proj_cs
        self.gt = tarProj.gt
        
    def GetGeometry(self): 
        #Get the necessary band information
        self.lins = self.datasource.RasterYSize
        self.cols = self.datasource.RasterXSize
        self.cellnull = self.layer.GetNoDataValue()
        self.celltype = gdal.GetDataTypeName(self.layer.DataType)
        self.projection = self.datasource.GetProjection()
        self.geotrans = self.datasource.GetGeoTransform()

        #Get the extent of the image
        self.ext=[]
        xarr=[0,self.cols]
        yarr=[0,self.lins]
        for px in xarr:
            for py in yarr:
                x=self.gt[0]+(px*self.gt[1])+(py*self.gt[2])
                y=self.gt[3]+(px*self.gt[4])+(py*self.gt[5])
                self.ext.append([x,y])
            yarr.reverse()
        self.bounds = (self.ext[0][0], self.ext[2][1], self.ext[2][0],self.ext[0][1])
        #Get the spatial resolution
        cellsize = [(self.ext[2][0]-self.ext[0][0])/self.cols, (self.ext[0][1]-self.ext[2][1])/self.lins] 
        if cellsize[0] != cellsize[1]:
            pass
        self.cellsize = cellsize[0]
 
    def ReadBand(self):
        self.NPBAND = self.layer.ReadAsArray()
        
        
    def ReadBlock(self, col,row,ncols,nrows):
        
        NPBLOCK = self.layer.ReadAsArray(col, row, ncols, nrows)
        self.NPBLOCK = NPBLOCK.ravel()


    def Flatten2DtoArr(self,structcodeD):
        self.GetGeometry()
        structcode  = structcodeD[self.cellType]

        bandList = list(chain.from_iterable(self.NPBAND))   
        self.bandArray = arr.array(structcode)
        self.bandArray.fromlist(bandList)
  
    def ColorTable(self,colortable):
        print ('    Setting color table (0 if success):'),
        print ( self.layer.SetColorTable(colortable) )
        
    def CloseLayer(self):
        #close the layer
        self.layer = None
 
   
def ESRICreateDSLayer(tarShpFPN,spatialRef,geomtype,layerid,fieldDefL):
    '''Creates a standard ESRI datasource and sets the layer and layerDef
    '''
    #Create an instance of datasource for the target
    tarDS = VectorDataSource()
    
    #Create the datasource
    tarDS.CreateESRIdatasource(tarShpFPN)
    
    #create an instance of layer
    tarLayer = VectorLayer()
    
    #link the datasource to the layer
    tarLayer.SetDS(tarDS)
    
    tarLayer.SetCSSpatialRef(spatialRef) 
      
    tarLayer.CreateOgrLayer(geomtype, str(layerid))
    
    tarLayer.GetLayerFields()
    
    for fieldDef in fieldDefL:
        
        tarLayer.AddFieldDef(fieldDef)
        
    return tarDS, tarLayer
            
def CreateESRIPolygonPtL(tarShpFPN, fieldDefL, ptL, proj_cs, layerid):
    ''' Create ESRI Polygon vector from list of points
    '''
    
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, 'polygon', layerid, fieldDefL)
    
    #create the geometry
    geom = ShapelyPolyGeom(ptL)
    
    geom.ShapelyToOgrGeom()
    
    tarFeat = ogrFeature(tarLayer)
    
    tarFeat.CreateOgrFeature(geom.ogrGeom, fieldDefL)
    
    tarDS.CloseDS()
    
def CreateESRIPolygonGeom(tarShpFPN, fieldDefL, geom, proj_cs, layerid):
    '''
    '''
    
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, 'polygon', layerid, fieldDefL)
    
    #create the geometry
    
    geom.ShapelyToOgrGeom()
    
    tarFeat = ogrFeature(tarLayer)
    
    tarFeat.CreateOgrFeature(geom.ogrGeom, fieldDefL)
    
    tarDS.CloseDS()
    
def ESRIOpenGetLayer(srcShpFPN,mode='read'):
    '''Opens a standard ESRI datasource and gets the layer and layerDef
    '''
    
    #Create an instance of datasource for the source data
    
    srcDS = VectorDataSource()
    
    #open the source data for reading
    
    if mode[0] == 'r':
    
        srcDS.OpenESRIshapeRead(srcShpFPN)
    
    else:
    
        srcDS.OpenESRIshapeEdit(srcShpFPN)
    
    #create a layer instance
    
    srcLayer = VectorLayer()
    
    #set the datasource of the layer
    
    srcLayer.SetDS(srcDS)
    
    #get the layer and the fieldDefintions 
    
    srcLayer.GetLayer()
    
    #get the spatialref
    
    srcLayer._GetSpatialRef()
    
    return srcDS,srcLayer



def ExtractFeaturesToNewDS(srcShpFPN,tarShpFPN,fieldname,valueLL,fieldDefL = False, dictL = False, key = False):
    '''
    '''
    
    #Open ESRI DS for reading
    srcDS,srcLayer = ESRIOpenGetLayer(srcShpFPN,'read')
    
    srcFieldDefL = srcLayer.fieldDefL
    
    if not fieldDefL:
        
        tarFieldDefL = srcFieldDefL
        
    else:
        
        tarFieldDefL = fieldDefL
        
    spatialRef = srcLayer.spatialRef

    #Create an instance of datasource for the target
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, spatialRef, srcLayer.geomtype, srcLayer.layerid, fieldDefL)
    
    #loop over all the features in the source file        
    for feature in srcLayer.layer:
        
        #create an instance of feature from the srcLayer   
        srcFeat = ogrFeature(srcLayer)
        
        #set the feature
        srcFeat.SetFeature(feature)
        
        #Set the fieldDef source from the source feature
        srcFeat.AddSourceValtoFieldDef( srcFieldDefL )
        
        extract = False
        
        for valueL in valueLL:
            
            for f in srcFieldDefL:
                
                if f.name == fieldname and f.source in valueL:
                    
                    extract = True
                    
                if extract:
                    
                    #create an instance of feature for the target
                    
                    tarFeat = ogrFeature(tarLayer) 
                    
                    tarFeat.CopyOgrFeature(feature,tarFieldDefL) 
                     
                    break #as the correct field was found  
                    
    srcDS.CloseDS()
    
    tarDS.CloseDS()
    
def GetFeatureBounds(srcShpFPN,fieldId):
    '''
    '''
    srcDS,srcLayer = ESRIOpenGetLayer(srcShpFPN,'read')
    
    boundsD = {}
    
    for feature in srcLayer.layer:
    
        srcFeat = ogrFeature(srcLayer)
        
        srcFeat.SetFeature(feature)
        
        fid = feature.GetField( str(fieldId) )
        
        geom = Geometry()
        
        geom.GeomFromFeature(feature)      
        
        # Always set the id with lower letters and hyphen to avoid mistakes
        boundsD[fid.lower().replace(' ','-')] = geom.shapelyGeom.bounds
    
    srcDS.CloseDS()
    
    return boundsD

def ReprojectBounds(bounds,cs_src,cs_tar):
    '''
    '''
    
    transform = osr.CoordinateTransformation(cs_src,cs_tar)
    
    ptL = []
    
    for pt in bounds:
    
        ptgeom = ShapelyPointGeom(pt)
        
        ptgeom.ShapelyToOgrGeom()
        
        #ptgeom.
        
        #point = ogr.CreateGeometryFromWkt("POINT (1120351.57 741921.42)")
        
        ptgeom.ogrGeom.Transform(transform)
        
        ptgeom.OgrGeomToShapely()
        
        ptcoord = list(ptgeom.shapelyGeom.coords)[0]
        
        ptL.extend([ptcoord[0],ptcoord[1]])     
    
    paramL = ['ullon','ullat','urlon','urlat','lrlon','lrlat','lllon','lllat']
    
    return dict(zip(paramL,ptL))

def GetVectorProjection(srcShpFPN):
    ''' Get the datasource and layer from a shape file
    '''
    
    filetype = os.path.splitext(srcShpFPN)[1]
    
    if filetype.lower() == '.shp':
    
        srcDS,srcLayer = ESRIOpenGetLayer(srcShpFPN,'read')
        
        srcLayer._GetSpatialRef()
        
        spatialRef = MjProj()
        
        if srcLayer.spatialRef != None:
        
            spatialRef.SetProj(srcLayer.spatialRef)
            
            spatialRef.ReadSpatialRef()
        
        else:
        
            spatialRef.epsg = False
        
        srcDS.CloseDS()
        
        return spatialRef
    
    else:
        
        errorstr = 'urecognized format for retrieveing spatial reference: %s' %(filetype)
        
        sys.exit(errorstr)
            
def RasterOpenGetFirstLayer(srcRastFPN,modeD):
    '''Opens a standard GDAL Raster datasource and gets the first layer
    '''
    #Create an instance of datasource for the source data
    srcDS = RasterDataSource()
    #open the source data for reading    
    if modeD['mode'] == 'read':
        srcDS.OpenGDALRead(srcRastFPN)
        #create a layer instance
        srcLayer = RasterLayer()
        #set the datasource of the layer
        srcLayer.SetDS(srcDS)
    elif modeD['mode'] == 'edit':
        srcDS.OpenGDALEdit(srcRastFPN)
        srcLayer = RasterLayer()
        srcLayer.SetDS(srcDS)
    else:
        print ( 'Can not understand edit mode', modeD['mode'] )
    srcLayer.GetLayer(1)
    #get the spatialref
    srcLayer.GetSpatialRef()
    srcLayer.GetGeometry()
    return srcDS,srcLayer 

def RasterCreateWithFirstLayer(dstRastFPN,layer):
    '''
    '''
    dstDS = RasterDataSource()
    
    dstDS.CreateGDALDS(dstRastFPN,layer)
    
    return dstDS
        
def GetRasterMetaData(srcRasterFPN): 
    '''
    '''
    
    modeD = {'mode':'read'}
    
    srcDS,srcLayer = RasterOpenGetFirstLayer(srcRasterFPN,modeD)                      
    
    #srcRast = RasterDataSource()
    spatialRef = MjProj()
    
    #TGTODO HERE I ASSUME THAT RASTER PROJ IS IN WKT FORMAT
    spatialRef.SetFromWKT(srcLayer.spatialRef)
    
    spatialRef.SetProj4()
    
    #spatialRef.SetProj(srcLayer.spatialRef)
    spatialRef.ReadSpatialRef()
    
    srcLayer.GetGeometry()
    
    #close the layer
    srcLayer.CloseLayer()
    srcDS.CloseDS()
    return spatialRef, srcLayer

def GetFeatureAttributeList(srcShpFPN, fieldL, idfield):
    '''
    '''

    srcDS, srcLayer = ESRIOpenGetLayer(srcShpFPN,'edit')
    
    fieldNameL = []
    
    fieldD = {}
        
    for srcFieldDef in srcLayer.fieldDefL:
    
        fieldNameL.append(srcFieldDef.name)
    
    for f in fieldL:
    
        if not f in fieldNameL:
        
            printstr = 'the requested field %s does not exists in the file %s' %(f, srcShpFPN)
            
            print ( printstr )
            
            return
    
    for feature in srcLayer.layer:
    
        values_list = [feature.GetField(str(j)) for j in fieldL]
        
        fid = feature.GetField(str(idfield))
        
        fieldD[fid] = dict(zip(fieldL,values_list))
    
    srcDS.CloseDS()
    
    return fieldD