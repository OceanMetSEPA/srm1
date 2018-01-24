# -*- coding: utf-8 -*-
"""
SRM1 Road Network

Created on Tue Oct 04 12:15:08 2016

@author: edward.barratt
"""
from os import path

import pandas as pd
import geopandas as gpd
from numpy import zeros_like, ones_like, reshape, array, append
from shapely.geometry import Point
#from VehicleCounts import VehicleCounts
#from EmissionFactors import EmissionFactorClass



class RoadNetwork(object):
  """
  Attributes:
    geoDataFrame
    emissionFactorName
    emissionFactors
    averageWindSpeed

  """
  __DispersionCoefficientsAll = {'NarrowCanyon':   {'A': 0.000488, 'B': -0.0308, 'C': 0.59, 'Alpha': 0},
                                 'WideCanyon':     {'A': 0.000325, 'B': -0.0205, 'C': 0.39, 'Alpha': 0.856},
                                 'OneSidedCanyon': {'A': 0.000500, 'B': -0.0316, 'C': 0.57, 'Alpha': 0},
                                 'NotACanyon':     {'A': 0.000310, 'B': -0.0182, 'C': 0.33, 'Alpha': 0.799}}
  __ImpactDistances = {'NarrowCanyon': 30,
                       'WideCanyon': 60,
                       'OneSidedCanyon': 30,
                       'NotACanyon': 60}
  __constants = {'B': 0.6, 'K': 100}
  __allowedTreeFactors = [1, 1.25, 1.5]

  def __init__(self, roadnetwork, averageWindSpeed=5,
               backgroundConcentration='AllZero',
               treefactorFieldName=None,
               NOxEmissionFieldName=None,
               NO2EmissionFieldName=None,
               otherChemsEmissionFieldNames=[None],
               isCanyonFieldName=None,
               roadWidthFieldName=None,
               canyonDepthFieldName=None):
    """

    required: isCanyonFieldName, roadWidthFieldName

    """

    self.__beingCreated = True

    # Get the defaults.
    self.treeFactorFieldName = treefactorFieldName
    self.NOxEmissionFieldName = NOxEmissionFieldName
    self.NO2EmissionFieldName = NO2EmissionFieldName
    self.otherChemsEmissionFieldNames = otherChemsEmissionFieldNames
    self.isCanyonFieldName = isCanyonFieldName
    self.roadWidthFieldName = roadWidthFieldName
    self.canyonDepthFieldName = canyonDepthFieldName

    self.averageWindSpeed = averageWindSpeed


    # Initiate a pandas dataframe with the correct columns.
    self.geoDataFrame = roadnetwork
    self.backgroundConcentration = backgroundConcentration
    self.__beingCreated = False
    #self.CalculateEmissions()      # Will create a column 'emissions_Pol' for each pollutant Pol.
    self.CalculateRoadConcentrations()  # Will create a column 'concentration_Pol' for each pollutant Pol.
    #self.bounds = self.GetBounds()

  @property
  def geoDataFrame(self):
    return self.__geoDataFrame

  @geoDataFrame.setter
  def geoDataFrame(self, v):
    print('geoDataFrame.setter called.')

    if isinstance(v, str):
      # Assume it's a shape file.
      v = gpd.read_file(v)
    elif isinstance(v, gpd.geodataframe.GeoDataFrame):
      pass
    else:
      raise ValueError('Unknown type for geoDataFrame.')

    # Some value checkers. Should probably add more.
    colNames = list(v)
    # Tree Factor.
    self.treeFactorFieldName, v = getFieldName(self.treeFactorFieldName, v, possibilities=['treefactor', 'tree'], defaultName='treefactor', defaultValue=1, name='tree factor')
    # Test values for tree factor
    Test = v[self.treeFactorFieldName]
    for Tu in Test.unique():
      if Tu not in self.__allowedTreeFactors:
        raise ValueError('Tree factor "{}" is not allowed.'.format(str(Tu)))

    # road Width
    self.roadWidthFieldName, v = getFieldName(self.roadWidthFieldName, v, possibilities=['roadwidth', 'width'], name='road width', required=True)
    # road Depth
    self.canyonDepthFieldName, v = getFieldName(self.canyonDepthFieldName, v, possibilities=['canyon', 'depth'], name='canyon depth', required=False)
    # is canyon
    self.isCanyonFieldName, v = getFieldName(self.isCanyonFieldName, v, possibilities=['is_canyon', 'iscanyon', 'roadclass'], name='canyon class', required=True)
    v = canyonDeciderAll(v, 'roadclass_c', self.isCanyonFieldName, self.roadWidthFieldName, self.canyonDepthFieldName)
    self.roadClassFieldNames = 'roadclass_c'

    # NOx Emissions
    avPols = {}
    if self.NOxEmissionFieldName is None:
      # None defined, test a few obvious possibilities.
      tests = ['nox (g/km)', 'nox emissions (g/km)', 'eft nox (g']
      gotNOx = False
      for colName in colNames:
        if colName.lower()[:10] in tests:
          print('NOx Emissions field name discovered: {}.'.format(colName))
          self.NOxEmissionFieldName = colName
          avPols['NOx'] = colName
          gotNOx = True
          break
      if not gotNOx:
        print('No NOx emissions found.')
    # NO2 Emissions
    if self.NO2EmissionFieldName is None:
      # None defined, test a few obvious possibilities.
      tests = ['no2 (g/km)', 'no2 emissi', 'eft no2 (g']
      gotNO2 = False
      for colName in colNames:
        if colName.lower()[:10] in tests:
          print('NO2 Emissions field name discovered: {}.'.format(colName))
          self.NO2EmissionFieldName = colName
          avPols['NO2'] = colName
          gotNO2 = True
          break
      if not gotNO2:
        print('No NO2 emissions found.')
      else:
        if not gotNOx:
          raise ValueError(("NO2 emissions detected, but no NOx. SRM1 cannot "
                            "calculate NO2 concentrations without both NOx "
                            "and NO2 emissions."))
    # other Chemicals
    if self.otherChemsEmissionFieldNames[0] is None:
      # None defined, test a few obvious possibilities.
      otherChemTests = {'PM10': ['pm10 (g/km', 'pm10 emiss', 'eft pm10 ('],
                        'PM25': ['pm25 (g/km', 'pm25 emiss', 'pm2.5 (g/k', 'pm2.5 emis', 'eft pm25 (', 'eft pm25 (', 'eft pm2.5']}
      gotOther = False
      self.otherChemsEmissionFieldNames = []
      for colName in colNames:
        for pol, tests in otherChemTests.items():
          if colName.lower()[:10] in tests:
            print('{} emissions field name discovered: {}.'.format(pol, colName))
            self.otherChemsEmissionFieldNames.append(colName)
            avPols[pol] = colName
            gotOther = True
      if (not gotNOx) and (not gotOther):
        raise ValueError("No chemical emission fields are recognised.")

    self.__geoDataFrame = v
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

  """
  def set_value(self, indices, col, values):
    changeableColumnsNoEffect = ['roadID', 'roadName']
    changeableColumnsReCalculate = ['roadClass', 'width', 'treeFactor',
                                    'speedClass', 'speed', 'stagnation']
    if col not in changeableColumnsNoEffect + changeableColumnsReCalculate:
      raise ValueError('You cannot change column "{}".'.format(col))
    else:
      if not isinstance(indices, list):
        indices = [indices]
      if not isinstance(values, list):
        values = [values]*len(indices)
      # Change the values
      for ii, ix in enumerate(indices):
        self.geoDataFrame.set_value(ix, col, values[ii])
      if col in changeableColumnsReCalculate:
        self.CalculateEmissions()
        self.CalculateRoadConcentrations()
   """

  def CalculatePointConcentrations(self, points, saveloc=None, head=False):

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


    if saveloc is not None:
      points.to_file(saveloc, driver='ESRI Shapefile', crs_wkt=crs_wkt)
    return points
    #points.apply(lambda row: self.CalculatePointConcentration(row), axis=1)


  def CalculatePointConcentration(self, point, pols=None, roads_impacting=None, closest_only=False):
    if pols is None:
      pols = self.availablePollutants

    if roads_impacting is None:
      pt, roads_impacting = self.RoadsImpactingPoint(point, closest_only=closest_only)

      DispFactor = (roads_impacting['dispersionCoefficient_C'] +
                    roads_impacting['dispersionCoefficient_B'] * roads_impacting['distance_point'] +
                    roads_impacting['dispersionCoefficient_A'] * roads_impacting['distance_point']**2)
      DispFactorB = roads_impacting['dispersionCoefficient_Alpha'] * roads_impacting['distance_point']**(-0.747)
      DispFactorB[roads_impacting['dispersionCoefficient_Alpha'] == 0] = 0
      DispFactor[roads_impacting['distance_point'] > 30] = DispFactorB[roads_impacting['distance_point'] > 30]
      DispFactor[roads_impacting['distance_point'] > roads_impacting['impactDistance']] = 0
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
    roads_int = roads_int[roads_int['distance_point'] < roads_int['impactDistance']].copy()
    roads_int = roads_int.sort_values(by=['distance_point'])

    if closest_only:
      roads_int = roads_int.head(1)

    pt['numRoads'] = len(roads_int.index)
    pt['Roads'] = ', '.join(list([str(x) for x in roads_int.index]))

    return pt, roads_int


  def CalculateRoadConcentrations(self, distance='RoadEdge', distanceAdd=0, distanceMultiply=1):
    if self.__beingCreated:
      return

    ImpactDistances = [self.__ImpactDistances[s] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['impactDistance'] = ImpactDistances
    DispersionCoefficientA = [self.__DispersionCoefficientsAll[s]['A'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispersionCoefficient_A'] = DispersionCoefficientA
    DispersionCoefficientB = [self.__DispersionCoefficientsAll[s]['B'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispersionCoefficient_B'] = DispersionCoefficientB
    DispersionCoefficientC = [self.__DispersionCoefficientsAll[s]['C'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispersionCoefficient_C'] = DispersionCoefficientC
    DispersionCoefficientAlpha = [self.__DispersionCoefficientsAll[s]['Alpha'] for s in list(self.geoDataFrame[self.roadClassFieldNames])]
    self.geoDataFrame['dispersionCoefficient_Alpha'] = DispersionCoefficientAlpha

    DD = self.geoDataFrame[self.roadWidthFieldName]/2
    if distance == 'RoadEdge':
      distance = DD
    elif distance == 'RoadCentre':
      distance = zeros_like(DD)
    elif len(distance) == 1:
      distance = distance * ones_like(DD)
    distance = distance * distanceMultiply + distanceAdd

    DispFactor = (self.geoDataFrame['dispersionCoefficient_C'] +
                  self.geoDataFrame['dispersionCoefficient_B'] * distance +
                  self.geoDataFrame['dispersionCoefficient_A'] * distance**2)
    DispFactorB = self.geoDataFrame['dispersionCoefficient_Alpha'] * distance**(-0.747)
    DispFactorB[self.geoDataFrame['dispersionCoefficient_Alpha'] == 0] = 0
    DispFactor[distance > 30] = DispFactorB[distance > 30]
    DispFactor[distance > self.geoDataFrame['impactDistance']] = 0
    #print distance
    DispFactor = DispFactor * 0.62 * 5 * self.geoDataFrame[self.treeFactorFieldName] / self.averageWindSpeed

    AvPols = self.availablePollutants

    # Do NOx first.
    if 'NOx' in AvPols.keys():
      self.geoDataFrame['concentration_NOx'] = (1000/(24*3600))*DispFactor*self.geoDataFrame[AvPols['NOx']] + self.backgroundConcentration['NOx']
    # then NO2
    if 'NO2' in AvPols.keys():
      fNO2 = self.geoDataFrame[AvPols['NO2']]/self.geoDataFrame[AvPols['NOx']]

      q = self.geoDataFrame['concentration_NOx']*(1-fNO2)
      w = q/(q+self.__constants['K'])
      e = self.__constants['B']*self.backgroundConcentration['O3']*w
      self.geoDataFrame['concentration_NOx'] = fNO2*self.geoDataFrame['concentration_NOx'] + e

    for Pol, FName in AvPols.items():
      if Pol in ['NOx', 'NO2']:
        continue
      self.geoDataFrame['concentration_{}'.format(Pol)] = (1000/(24*3600))*DispFactor*self.geoDataFrame[FName] + self.backgroundConcentration[Pol]

def getFieldName(initialFN, data, possibilities=[], defaultName=None, defaultValue=None, required=False, name='unnamed'):
  colNames = list(data)
  if initialFN is None:
    got = False
    for colName in colNames:
      if colName.lower() in possibilities:
        print('Field name discovered for {}: {}.'.format(name, colName))
        FName = colName
        got = True
        break
    if not got:
      if defaultName is None:
        if required:
          raise ValueError('No field name assigned for {}'.format(name))
        else:
          FName = None
      else:
        print('No field name assigned for {}. Creating field {} with initial value {}.'.format(name, defaultName, defaultValue))
        if defaultName in colNames:
          raise ValueError('A field named {} already exists.'.format(defaultName))
        else:
          data[defaultName] = defaultValue
          FName = defaultName
  else:
    if initialFN not in colNames:
      raise ValueError("The assigned field name for {}, '{}', is not present in the data.".format(name, initialFN))
    else:
      FName = initialFN
  return FName, data


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
  if IsCanyon.lower() in ['n', 'no', 'notacanyon', 'not a canyon']:
    return 'NotACanyon'
  elif IsCanyon.lower() in ['1-sided', '1sided', 'onesidedcanyon', 'one sided canyon']:
    return 'OneSidedCanyon'
  elif IsCanyon.lower() in ['widecanyon', 'wide canyon']:
    return 'WideCanyon'
  elif IsCanyon.lower() in ['narrowcanyon', 'narrow canyon']:
    return 'NarrowCanyon'
  elif IsCanyon.lower() in ['y', 'yes']:
    if all([CanyonDepth, CanyonWidth]):
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
    raise ValueError('IsCanyon value {} is not understood.'.format(IsCanyon))


if __name__ == '__main__':
  rdsshpfile = ("\\\sepa-fp-01\DIR SCIENCE\EQ\Oceanmet\Projects\\air\CAFS\\"
             "Glasgow\Eddy\SourceFiles\GlasgowTracsisCentreT_wEmissions2017.shp")
  background = {'NO2': 25.5, 'NOx': 42.0, 'O3':38.6, 'PM10': 0, 'PM25': 0}

  RN = RoadNetwork(rdsshpfile)
  RN.backgroundConcentration = background

  ptsshpfile = ("\\\sepa-fp-01\DIR SCIENCE\EQ\Oceanmet\Projects\\air\CAFS\\"
                "Glasgow\Eddy\SourceFiles\outputPoints.shp")


  CP = RN.CalculatePointConcentrations(ptsshpfile)
  print(CP)
