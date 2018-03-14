# -*- coding: utf-8 -*-
"""
SRM1 Road Network

Created on Tue Oct 04 12:15:08 2016
@author: edward.barratt
"""
#from os import path

import pandas as pd
import geopandas as gpd
from numpy import zeros_like, ones_like
from shapely.geometry import Point
#from VehicleCounts import VehicleCounts
#from EmissionFactors import EmissionFactorClass

__notACanyonOptions = ['n', 'no', 'notacanyon', 'not a canyon', 'no canyon']
__oneSidedCanyonOptions = ['1-sided', '1sided', 'onesidedcanyon', 'one sided canyon']
__wideCanyonOptions = ['widecanyon', 'wide canyon']
__narrowCanyonOptions = ['narrowcanyon', 'narrow canyon']
__decideCanyonOptions = ['y', 'yes', '2-sided', 'complex canyon']

class RoadNetwork(object):
  """
  An SRM1 RoadNetwork object.

  Args:
    geoDataFrame(str, or geopandas dataframe): A geopandas dataframe containing
                            the details of the road network. It must contain a
                            'geometry' field, and a number of other fields
                            detailed below. If provided as a string, then the
                            string will be assumed to be a path to a shapefile
                            and the geoDataFrame will be read from the shape file.
    treeFactorFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the tree factor for the
                            road segment. This field must only contain values 1,
                            1.25, or 1.5. If not provided then the user will be
                            invited to specify a field, or to create one with
                            default value of 1.
    isCanyonFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents whether or not the road
                            segment is a canyon. This column can contain any of
                            the following values:  ['n', 'no', 'notacanyon',
                            'not a canyon', '1-sided', '1sided', 'onesidedcanyon',
                            'one sided canyon', 'widecanyon', 'wide canyon',
                            'narrowcanyon', 'narrow canyon', 'y', 'yes']. If the
                            value is either 'y' or 'yes' then the nature of the
                            canyon (wide, narrow, or not a canyon) will be
                            determined by the roadWidth field and the
                            canyonDepth field. If not provided then the user
                            will be invited to specify a field.
    roadWidthFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the road width, in
                            meters. The field should contain numeric values.
                            These values are used to calculate road side
                            concentrations and to determine the canyon type of
                            roads (isCanyonFieldName). If not provided then the
                            user will be invited to specify a field.
    canyonDepthFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the canyon depth, in
                            meters. The field should contain numeric values.
                            These values are used to determine the canyon type
                            of roads (isCanyonFieldName). If not provided then the
                            user will be invited to specify a field, or to
                            create one with default value of 1; this last option
                            is not recomended if the isCanyon field contains
                            values 'y' or 'yes', see above for details.
    NOxEmissionFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the NOx emissions,
                            in g/km. The field should contain numeric values. If
                            not provided then the user will be invited to
                            specify a field.
    NO2EmissionFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the NO2 emissions,
                            in g/km. The field should contain numeric values. If
                            not provided then the user will be invited to
                            specify a field. Note that both NO2 and NOx
                            emissions must be specified in order to calculate
                            NO2 concentrations.
    PM10EmissionFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the PM10 emissions,
                            in g/km. The field should contain numeric values. If
                            not provided then the user will be invited to
                            specify a field.
    PM25EmissionFieldName(str, optional): The name of the field within the
                            geoDataFrame that represents the PM2.5 emissions,
                            in g/km. The field should contain numeric values. If
                            not provided then the user will be invited to
                            specify a field.
    backgroundConcentrations(dict, or the string 'AllZero', optional): A
                            dictionary containing the background concentrations
                            of the pollutants of interest. If NO2 concentrations
                            are to be calculated then a background O3
                            concentration is also required. Units are ug/m3. If
                            'AllZero', then all background concentrations will
                            be zero. Default 'AllZero'.
  """

  # Define some internal defaults for the class
  __DispersionCoefficientsAll = {'NarrowCanyon':   {'A': 0.000488, 'B': -0.0308, 'C': 0.59, 'Alpha': 0},
                                 'WideCanyon':     {'A': 0.000325, 'B': -0.0205, 'C': 0.39, 'Alpha': 0.856},
                                 'OneSidedCanyon': {'A': 0.000500, 'B': -0.0316, 'C': 0.57, 'Alpha': 0},
                                 'NotACanyon':     {'A': 0.000310, 'B': -0.0182, 'C': 0.33, 'Alpha': 0.799}}
  __impactDists = {'NarrowCanyon': 30,
                   'WideCanyon': 60,
                   'OneSidedCanyon': 30,
                   'NotACanyon': 60}
  __constants = {'B': 0.6, 'K': 100}
  __allowedTreeFactors = [1, 1.25, 1.5]
  __autoCreatedFields = {}
  __maxRoadsPerCalcPoint = 2

  def __init__(self, roadnetwork,
               averageWindSpeed=5,
               backgroundConcentration='AllZero',
               treeFactorFieldName=None,
               NOxEmissionFieldName=None,
               NO2EmissionFieldName=None,
               PM10EmissionFieldName=None,
               PM25EmissionFieldName=None,
               isCanyonFieldName=None,
               roadWidthFieldName=None,
               canyonDepthFieldName=None):

    self.__beingCreated = True

    # Get the defaults.
    self.treeFactorFieldName = treeFactorFieldName
    self.NOxEmissionFieldName = NOxEmissionFieldName
    self.NO2EmissionFieldName = NO2EmissionFieldName
    self.PM10EmissionFieldName = PM10EmissionFieldName
    self.PM25EmissionFieldName = PM25EmissionFieldName
    self.isCanyonFieldName = isCanyonFieldName
    self.roadWidthFieldName = roadWidthFieldName
    self.canyonDepthFieldName = canyonDepthFieldName

    self.averageWindSpeed = averageWindSpeed


    # Initiate a pandas dataframe with the correct columns.
    self.geoDataFrame = roadnetwork
    self.backgroundConcentration = backgroundConcentration
    self.__beingCreated = False
    #self.CalculateEmissions()      # Will create a column 'emissions_Pol' for each pollutant Pol.
    self.CalculateRoadConcentrations()  # Will create a column 'conc_Pol' for each pollutant Pol.
    #self.bounds = self.GetBounds()

  @property
  def geoDataFrame(self):
    return self.__geoDataFrame

  @geoDataFrame.setter
  def geoDataFrame(self, v):
    if not self.__beingCreated:
      raise IOError('geoDataFrame cannot be reset.')
    if isinstance(v, str):
      # Assume it's a shape file.

      prj_file = v.replace('.shp', '.prj')
      crs_wkt = [l.strip() for l in open(prj_file,'r')][0]
      self.crs_wkt = crs_wkt
      v = gpd.read_file(v)
    elif isinstance(v, gpd.geodataframe.GeoDataFrame):
      pass
    else:
      print(v)
      print(type(v))
      raise ValueError('Unknown type for geoDataFrame.')

    # make sure that required fields are recognised.
    colNames = list(v)
    self.originalColNames = colNames.copy()
    # Tree Factor.
    self.treeFactorFieldName, v = getFieldName(v, initialFN=self.treeFactorFieldName, options=colNames, defaultName='treefactor', defaultValue=1, name='tree factor', canIgnore=False)
    if self.treeFactorFieldName in colNames:
      colNames.remove(self.treeFactorFieldName)
    else:
      self.__autoCreatedFields['treefactor'] = 1

    # is canyon
    self.isCanyonFieldName, v = getFieldName(v, initialFN=self.isCanyonFieldName, options=colNames, name='canyon class', defaultName='canyon class', defaultValue='narrow canyon', canIgnore=False)
    if self.isCanyonFieldName in colNames:
      colNames.remove(self.isCanyonFieldName)
    # road Width
    self.roadWidthFieldName, v = getFieldName(v, initialFN=self.roadWidthFieldName, options=colNames, name='road width', canIgnore=False)
    if self.roadWidthFieldName in colNames:
      colNames.remove(self.roadWidthFieldName)
    # road Depth
    self.canyonDepthFieldName, v = getFieldName(v, initialFN=self.canyonDepthFieldName, options=colNames, name='canyon depth', canIgnore=True)
    if self.canyonDepthFieldName in colNames:
      colNames.remove(self.canyonDepthFieldName)
    v = canyonDeciderAll(v, 'roadclass_c', self.isCanyonFieldName, self.roadWidthFieldName, self.canyonDepthFieldName)
    self.roadClassFieldNames = 'roadclass_c'

    avPols = {}
    # NOx Emissions
    self.NOxEmissionFieldName, v = getFieldName(v, initialFN=self.NOxEmissionFieldName, options=colNames, name='NOx Emission Rate (g/km)', canIgnore=True)
    if self.NOxEmissionFieldName in colNames:
      colNames.remove(self.NOxEmissionFieldName)
      avPols['NOx'] = self.NOxEmissionFieldName
      # NO2 Emissions, can't do without NOx
      self.NO2EmissionFieldName, v = getFieldName(v, initialFN=self.NO2EmissionFieldName, options=colNames, name='NO2 Emission Rate (g/km)', canIgnore=True)
      if self.NO2EmissionFieldName is not None:
        avPols['NO2'] = self.NO2EmissionFieldName
    # PM10 Emissions
    self.PM10EmissionFieldName, v = getFieldName(v, initialFN=self.PM10EmissionFieldName, options=colNames, name='PM10 Emission Rate (g/km)', canIgnore=True)
    if self.PM10EmissionFieldName in colNames:
      colNames.remove(self.PM10EmissionFieldName)
      avPols['PM10'] = self.PM10EmissionFieldName
   # PM25 Emissions
    self.PM25EmissionFieldName, v = getFieldName(v, initialFN=self.PM25EmissionFieldName, options=colNames, name='PM2.5 Emission Rate (g/km)', canIgnore=True)
    if self.PM25EmissionFieldName in colNames:
      colNames.remove(self.PM25EmissionFieldName)
      avPols['PM25'] = self.PM25EmissionFieldName

    if len(avPols.keys()) == 0:
      raise ValueError("No chemical emission fields are recognised.")

    self.checkValues(v)
    self.__geoDataFrame = v
    print('{} roads imported.'.format(len(v.index)))
    self.availablePollutants = avPols
    #self.bounds = self.geoDataFrame.total_bounds

  @property
  def averageWindSpeed(self):
    return self.__averageWindSpeed

  @averageWindSpeed.setter
  def averageWindSpeed(self, v):
    self.__averageWindSpeed = v
    self.CalculateRoadConcentrations()

  @property
  def backgroundConcentration(self):
    return self.__backgroundConcentration

  @backgroundConcentration.setter
  def backgroundConcentration(self, v):
    if v == 'AllZero':
      v = {}
      for Pol in self.availablePollutants.keys():
        v[Pol] = 0
      v['O3'] = 0
    self.__backgroundConcentration = v
    self.CalculateRoadConcentrations()

  ## Some helper functions for setting up the geodataframe.

  def checkValues(self, geodataframe):
    """
    Check the values within key columns of a geodataframe.
    """
    # Test values for tree factor
    Test = geodataframe[self.treeFactorFieldName]
    for Tu in Test.unique():
      if Tu not in self.__allowedTreeFactors:
        raise ValueError('Tree factor "{}" is not allowed.'.format(str(Tu)))
    # Should add more tests.

  def CalculatePointConcentrations(self, points, saveloc=None, head=False):
    """
    Calculate the pollutant concentrations at a set of points.

    Args:
      points(str, or geopandas dataframe): A geopandas dataframe containing
                            the points where the concentration is to be
                            calculated. It must contain a 'geometry' field.
                            If provided as a string, then the string will be
                            assumed to be a path to a shapefile and the
                            geoDataFrame will be read from the shape file.
      saveloc(str, optional): The path where the output shapefile should be
                            saved. The output shape file will have the same
                            fields as the input geodataframe, plus a 'PEC_XX'
                            and 'PC_XX' field for each pollutant. If not
                            specified the saveloc will be the same as the input
                            shapefile with '_conc' appended to the filename, or
                            if points was a geopandas dataframe and not a
                            shapefile, then the points will not be saved.

    Outputs:

    """
    if isinstance(points, str):
      # Assume it's a shape file.
      prj_file = points.replace('.shp', '.prj')
      if saveloc is None:
        saveloc = points.replace('.shp', '_conc.shp')
      crs_wkt = [l.strip() for l in open(prj_file,'r')][0]
      points = gpd.read_file(points)
    elif isinstance(points, gpd.geodataframe.GeoDataFrame):
      crs_wkt = None
    else:
      raise ValueError('Unknown type for geoDataFrame.')

    if head:
      points = points.head()

    pt_first = self.CalculatePointConcentration([points.iloc[0]['x'], points.iloc[0]['y']])
    colnames = list(pt_first)
    colnames.remove('geometry')
    colnames.remove('Roads')

    for colname in colnames:
      points[colname] = 0
    points['Roads'] = '-'

    for ri, row in points.iterrows():
      print('Calculating concentration at point {} of {}.'.format(ri, len(points.index)), end='\r', flush=True)
      pt_conc = self.CalculatePointConcentration([row['x'], row['y']])
      for colname in colnames:
        points.iloc[ri, points.columns.get_loc(colname)] = pt_conc[colname].values
      points.iloc[ri, points.columns.get_loc('Roads')] = ', '.join(pt_conc['Roads'].values)
    print('Calculation of point concentrations complete.')

    # Add the background concentration
    for pol in self.availablePollutants.keys():
      points['PEC_{}'.format(pol)] = points['PC_{}'.format(pol)] + self.backgroundConcentration[pol]

    if saveloc is not None:
      points.to_file(saveloc, driver='ESRI Shapefile', crs_wkt=crs_wkt)
    return points
    #points.apply(lambda row: self.CalculatePointConcentration(row), axis=1)


  def CalculatePointConcentration(self, point, pols=None, roads_impacting=None, closest_only=False):
    if pols is None:
      pols = self.availablePollutants

    if roads_impacting is None:
      pt, roads_impacting = self.RoadsImpactingPoint(point, closest_only=closest_only)

      DispFactor = (roads_impacting['dispC_C'] +
                    roads_impacting['dispC_B'] * roads_impacting['distance_point'] +
                    roads_impacting['dispC_A'] * roads_impacting['distance_point']**2)
      DispFactorB = roads_impacting['dispC_Ap'] * roads_impacting['distance_point']**(-0.747)
      DispFactorB[roads_impacting['dispC_Ap'] == 0] = 0
      DispFactor[roads_impacting['distance_point'] > 30] = DispFactorB[roads_impacting['distance_point'] > 30]
      DispFactor[roads_impacting['distance_point'] > roads_impacting['impactDist']] = 0
      DispFactor = DispFactor * 0.62 * 5 * roads_impacting[self.treeFactorFieldName] / self.averageWindSpeed
      roads_impacting['FullFactor'] = (1000/(24.0*3600))*DispFactor


    # Do NOx first.
    if 'NOx' in pols.keys():
      c_NOx = roads_impacting['FullFactor']*roads_impacting[pols['NOx']]
      pt['PC_NOx'] = c_NOx.sum()
    # then NO2
    if 'NO2' in pols.keys():
      fNO2 = roads_impacting[pols['NO2']]/roads_impacting[pols['NOx']]

      q = c_NOx*(1-fNO2)
      w = q/(q+self.__constants['K'])
      e = self.__constants['B']*self.backgroundConcentration['O3']*w
      v = fNO2*c_NOx + e
      pt['PC_NO2'] = v.sum()

    for pol, colname in pols.items():
      if pol in ['NOx', 'NO2']:
        continue
      else:
        PCName = 'PC_{}'.format(pol)
        v = roads_impacting['FullFactor']*roads_impacting[colname]
        pt[PCName] = v.sum()
    return pt

  def RoadsImpactingPoint(self, point, closest_only=False):
    """
    point should be a two item column, x and y, in the same coordinate system
    as the geoDataFrame.
    """

    # Get the point as a geodataframe.
    pt = pd.DataFrame([Point(point[0], point[1])], columns=['geometry'])
    pt = gpd.GeoDataFrame(pt, geometry='geometry')

    # Get the distance of all roads to this point.
    roads_int = self.geoDataFrame.copy()
    roads_int['distance_point'] = roads_int.geometry.apply(lambda g: pt.distance(g))
    roads_int = roads_int[roads_int['distance_point'] < roads_int['impactDist']].copy()
    roads_int = roads_int.sort_values(by=['distance_point'])

    if closest_only:
      roads_int = roads_int.head(1)
    else:
      # Keep only the first few roads (usually 2, controlled by __maxRoadsPerCalcPoint)
      # This aims to prevent double counting. At a cross roads it probably makes
      # sense to count two roads, but not all 4, for example. There will still be
      # some double counting where two paralel segments of the same road meet.
      roads_int = roads_int.head(self.__maxRoadsPerCalcPoint)

    pt['numRoads'] = len(roads_int.index)
    pt['Roads'] = ', '.join(list([str(x) for x in roads_int.index]))

    return pt, roads_int


  def CalculateRoadConcentrations(self, distance='RoadEdge', distanceAdd=0, distanceMultiply=1):
    if self.__beingCreated:
      return

    impactDists = [self.__impactDists[s] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['impactDist'] = impactDists
    DispersionCoefficientA = [self.__DispersionCoefficientsAll[s]['A'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispC_A'] = DispersionCoefficientA
    DispersionCoefficientB = [self.__DispersionCoefficientsAll[s]['B'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispC_B'] = DispersionCoefficientB
    DispersionCoefficientC = [self.__DispersionCoefficientsAll[s]['C'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispC_C'] = DispersionCoefficientC
    DispersionCoefficientAlpha = [self.__DispersionCoefficientsAll[s]['Alpha'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispC_Ap'] = DispersionCoefficientAlpha

    DD = self.geoDataFrame[self.roadWidthFieldName]/2
    if distance == 'RoadEdge':
      distance = DD
    elif distance == 'RoadCentre':
      distance = zeros_like(DD)
    elif len(distance) == 1:
      distance = distance * ones_like(DD)
    distance = distance * distanceMultiply + distanceAdd

    DispFactor = (self.geoDataFrame['dispC_C'] +
                  self.geoDataFrame['dispC_B'] * distance +
                  self.geoDataFrame['dispC_A'] * distance**2)
    DispFactorB = self.geoDataFrame['dispC_Ap'] * distance**(-0.747)
    DispFactorB[self.geoDataFrame['dispC_Ap'] == 0] = 0
    DispFactor[distance > 30] = DispFactorB[distance > 30]
    DispFactor[distance > self.geoDataFrame['impactDist']] = 0
    #print distance
    DispFactor = DispFactor * 0.62 * 5 * self.geoDataFrame[self.treeFactorFieldName] / self.averageWindSpeed

    AvPols = self.availablePollutants

    # Do NOx first.
    if 'NOx' in AvPols.keys():
      #self.geoDataFrame['conc_NOx']
      PCNOX = (1000/(24*3600))*DispFactor*self.geoDataFrame[AvPols['NOx']]
    # then NO2
    if 'NO2' in AvPols.keys():
      fNO2 = self.geoDataFrame[AvPols['NO2']]/self.geoDataFrame[AvPols['NOx']]

      q = PCNOX*(1-fNO2)
      w = q/(q+self.__constants['K'])
      e = self.__constants['B']*self.backgroundConcentration['O3']*w
      self.geoDataFrame['conc_NO2'] = fNO2*PCNOX + e + self.backgroundConcentration['NO2']
    if 'NOx' in AvPols.keys():
      self.geoDataFrame['conc_NOx'] = PCNOX + self.backgroundConcentration['NOx']
    for Pol, FName in AvPols.items():
      if Pol in ['NOx', 'NO2']:
        continue
      self.geoDataFrame['conc_{}'.format(Pol)] = (1000/(24*3600))*DispFactor*self.geoDataFrame[FName] + self.backgroundConcentration[Pol]

  def addRoads(self, v):

    if isinstance(v, str):
      # Assume it's a shape file.
      v = gpd.read_file(v)
    elif isinstance(v, gpd.geodataframe.GeoDataFrame):
      pass
    else:
      raise ValueError('Unknown type for geoDataFrame.')

    print('Adding new roads to existing road network.')
    # Check which columns match those columns from the existing network.
    newColNames = list(v)
    origColNames = self.originalColNames

    matches = []
    origmissing = []
    newmissing = []
    for cn in origColNames:
      if cn in newColNames:
        matches.append(cn)
      else:
        origmissing.append(cn)
    for cn in newColNames:
      if cn not in matches:
        newmissing.append(cn)

    doNotSkip = [self.treeFactorFieldName, self.NOxEmissionFieldName,
                 self.NO2EmissionFieldName, self.PM10EmissionFieldName,
                 self.PM25EmissionFieldName, self.isCanyonFieldName,
                 self.roadWidthFieldName, self.canyonDepthFieldName]

    print('The following columns exist in both road networks, and will be linked.')
    print(', '.join(matches))
    print('')
    for cn in origmissing:
      if cn in doNotSkip:
        msg = ('Column "{}" exists in the original road network but is missing '
               'in the new one.\nPlease specify which column, to match it '
               'with.').format(cn)
        fName, v = getFieldName(v, options=newmissing, canIgnore=False, questionString=msg)
      else:
        msg = ('Column "{}" exists in the original road network but is missing '
               'in the new one.\nPlease specify which column, if any, to match it '
               'with.').format(cn)
        fName, v = getFieldName(v, options=newmissing, canIgnore=True, questionString=msg)
      if fName in newmissing:
        newmissing.remove(fName)

      print('')
    if len(newmissing):
      print('The following columns exist only in the new road network, and will be ignored.')
      print(', '.join(newmissing))
      print('')

    # Add the columns that were automatically added by the geodataframe setter.
    for key, value in self.__autoCreatedFields.items():
      v[key] = value

    self.checkValues(v)

    vall = self.geoDataFrame.copy()
    vall = vall.append(v, ignore_index=True)
    vall = canyonDeciderAll(vall, 'roadclass_c', self.isCanyonFieldName, self.roadWidthFieldName, self.canyonDepthFieldName)
    self.__geoDataFrame = vall
    print('{} roads imported. There are now {} roads in total.'.format(len(v.index), len(vall.index)))
    self.CalculateRoadConcentrations()

  def toFile(self, saveloc):
    """
    Save the geodataframe to a shapefile.
    This will save the complete geodataframe, with the calculated concentration
    fields, etc.

    args:
      saveloc (str): The path where the shape file is to be saved.
    """
    self.geoDataFrame.to_file(saveloc, crs_wkt=self.crs_wkt)

## Other helper functions
def canyonDeciderAll(DF, newColName, IsCanyonFN, CanyonDepthFN, CanyonWidthFN):
  if all([CanyonDepthFN, CanyonWidthFN]):
    DF[newColName] = DF.apply(lambda row: canyonDecider(row[IsCanyonFN],
                                                        row[CanyonDepthFN],
                                                        row[CanyonWidthFN]),
                                          axis=1)
  else:
    DF[newColName] = DF.apply(lambda row: canyonDecider(row[IsCanyonFN], None, None), axis=1)
  return DF

def canyonDecider(IsCanyon, CanyonDepth, CanyonWidth):
  if IsCanyon.lower() in __notACanyonOptions:
    return 'NotACanyon'
  elif IsCanyon.lower() in __oneSidedCanyonOptions:
    return 'OneSidedCanyon'
  elif IsCanyon.lower() in __wideCanyonOptions:
    return 'WideCanyon'
  elif IsCanyon.lower() in __narrowCanyonOptions:
    return 'NarrowCanyon'
  elif IsCanyon.lower() in __decideCanyonOptions:
    if CanyonDepth < 1:
      return 'NotACanyon'
    elif all([CanyonDepth, CanyonWidth]):
      # If we have both canyon depth and canyon width.
      decider = CanyonWidth/(2*CanyonDepth)
      if decider >= 3:
        return 'NotACanyon'
      elif decider < 1.5:
        return 'NarrowCanyon'
      else: # Between 1.5 and 3
        return 'WideCanyon'
    else:
      # Default is NarrowCanyon
      return 'NarrowCanyon'
  else:
    raise ValueError('IsCanyon value "{}" is not understood.'.format(IsCanyon))

def getFieldName(data, initialFN=None, options=None, defaultName=None, defaultValue=None, name='unnamed', canIgnore=True, questionString=None):

  if options is None:
    options = list(data)
  if initialFN is not None:
    if initialFN in options:
      # The initial suggestion works. Nothing more to do.
      return initialFN, data

  if questionString is None:
    questionString = 'Please select a column to use for field "{}"'.format(name)
  print(questionString)

  # Need to ask for suggestions.
  if canIgnore:
    print('Choose from the following:')
    optionNums = [str(x) for x in range(len(options)+1)]
  else:
    print('This column cannot be ignored.')
    print('Choose from the following:')
    optionNums = [str(x+1) for x in range(len(options))]
  for ci, cn in enumerate(options):
    print('{:3d}: {:18s}'.format(ci+1, cn), end='')
    nlneeded = True
    if (ci+1)%3 == 0:
      print('')
      nlneeded = False
  if canIgnore:
    print('  0: Ignore this field')
    nlneeded = False

  defaultNum = -99
  if defaultName is not None:
    defaultNum = len(options)+1
    print('{:3d}: Create a field named "{}" with default value {}.'.format(defaultNum, defaultName, defaultValue))
    nlneeded = False
    optionNums.append(str(defaultNum))
  if nlneeded:
    print('')
  num = input("Match column number? ")
  while num not in optionNums:
    print('Option not recognised.')
    num = input("Match column number? ")
  num = int(num)

  if num == defaultNum:
    data[defaultName] = defaultValue
    fName = defaultName
  elif num == 0:
    fName = None
  else:
    fName = options[num-1]

  return fName, data

if __name__ == '__main__':
  # Glasgow
  #rdsshpfile = ("C:\\Users\\edward.barratt\\Documents\\Modelling\\CAFS\\Glasgow\\"
  #              "Glasgow_AllRoads_wEmissions2017.shp")
  #background = {'NO2': 25.5, 'NOx': 42.0, 'O3':38.6, 'PM10': 0, 'PM25': 0}
  #
  # Start with the Central Glasgow file
  #RN = RoadNetwork(rdsshpfile, isCanyonFieldName='IsCanyon',
  #                 roadWidthFieldName='WIDTH', canyonDepthFieldName='CANYON',
  #                 averageWindSpeed=4)
  #RN.backgroundConcentration = background
  #
  #ptsshpfile = ("\\\sepa-fp-01\DIR SCIENCE\EQ\Oceanmet\Projects\\air\CAFS\\"
  #              "Glasgow\Eddy\SourceFiles\outputPoints.shp")

  #RN.CalculateRoadConcentrations(distance='RoadCentre', distanceAdd=12, distanceMultiply=1)
  #
  #saveloc = ("\\\sepa-fp-01\DIR SCIENCE\EQ\Oceanmet\Projects\\air\CAFS\\"
  #           "Glasgow\Eddy\SourceFiles\GlasgowDutch2017.shp")

  # Edinburgh
  rdsshpfile = ("C:\\Users\\edward.barratt\\Documents\\Modelling\\CAFS\\Edinburgh\\"
                "Roads\\RoadsForEMIT\\171129\\171129_EdinburghRoads_ForADMS_v6_JoinTo171128_TrafficData_CW_wEmissions2017.shp")
  background = {'NO2': 20.0, 'NOx': 28.0, 'O3':46, 'PM10': 10, 'PM25': 7}

  # Start with the Central Glasgow file
  RN = RoadNetwork(rdsshpfile, isCanyonFieldName='TYPE_CANYO',
                   roadWidthFieldName='WIDTH', canyonDepthFieldName='CANYON',
                   averageWindSpeed=4)
  RN.backgroundConcentration = background

  ptsshpfile = ("C:\\Users\\edward.barratt\\Documents\\Modelling\\CAFS\\"
                "Edinburgh\\Points\\reducedPoints.shp")

  #RN.CalculateRoadConcentrations(distance='RoadCentre', distanceAdd=12, distanceMultiply=1)
  #
  saveloc = ("C:\\Users\\edward.barratt\\Documents\\Modelling\\CAFS\\Edinburgh\\"
                "Roads\\RoadsForEMIT\\EdinburghDutch2017.shp")


  RN.toFile(saveloc)
  CP = RN.CalculatePointConcentrations(ptsshpfile)

