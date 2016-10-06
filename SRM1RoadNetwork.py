# -*- coding: utf-8 -*-
"""
SRM1 Road Network

Created on Tue Oct 04 12:15:08 2016

@author: edward.barratt
"""
import geopandas as gpd
from numpy import zeros_like, ones_like, reshape, array, append
from VehicleCounts import VehicleCounts
from EmissionFactors import EmissionFactorClass

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
  __PotentialVehMixes = {'DFT1': ['Motorcycle', 'PedalCycle', 'CarsTaxis', 'BusOrCoach', 'LGV', 'RHGV2Ax', 'RHGV3Ax', 'RHGV4or5Ax', 'AHGV3or4Ax', 'AHGV5Ax', 'AHGV6omAx'],
                         'DFT2': ['MCycle', 'Car', 'Bus', 'LGV', 'RHGV_2X', 'RHGV_3X', 'RHGV_4X', 'AHGV_34X', 'AHGV_5X', 'AHGV_6X'],
                         'GGow': ['MC', 'PC', 'CAR', 'PSV', 'LGV', 'OGV1', 'OGV2']}
  __allowedTreeFactors = [1, 1.25, 1.5]
  __columnNames = ['roadID', 'roadName', 'roadClass', 'speedClass', 'speed',
                   'width', 'stagnation', 'treeFactor', 'geometry', 'vehicleCounts']

  def __init__(self, V, emissionFactorName='NAEI', emissionFactors='Get',
                     averageWindSpeed=5, backgroundConcentration='AllZero'):
    self.__beingCreated = True

    # Initiate a pandas dataframe with the correct columns.
    self.geoDataFrame = gpd.GeoDataFrame(V)

    self.emissionFactorName = emissionFactorName
    self.emissionFactors = emissionFactors
    self.averageWindSpeed = averageWindSpeed
    self.backgroundConcentration = backgroundConcentration

    self.__beingCreated = False
    self.CalculateEmissions()      # Will create a column 'emissions_Pol' for each pollutant Pol.
    self.CalculateConcentrations() # Will create a column 'concentration_Pol' for each pollutant Pol.
    self.bounds = self.GetBounds()




  @property
  def geoDataFrame(self):
    return self.__geoDataFrame

  @geoDataFrame.setter
  def geoDataFrame(self, v):
    print 'geoDataFrame.setter called'
    self.__geoDataFrame = v
    # Some value checkers. Should probably add more.
    Test = self.geoDataFrame['treeFactor']
    for Tu in Test.unique():
      if Tu not in self.__allowedTreeFactors:
        raise ValueError('Tree factor "{}" is not allowed.'.format(str(Tu)))
    self.bounds = self.GetBounds()

  @property
  def emissionFactors(self):
    return self.__emissionFactors

  @emissionFactors.setter
  def emissionFactors(self, v):
    if v == 'Get':
      # Get the emission factors for the specified emission factor name.
      EF = EmissionFactorClass(Name=self.emissionFactorName)
    elif isinstance(v, EmissionFactorClass):
      EF = v
    elif isinstance(v, dict):
      EF = EmissionFactorClass(Factors=v)
    else:
      raise ValueError("emissionFactors should be an instance of emissionFactors.EmissionFactorClass.")
    self.__emissionFactors = EF
    self.CalculateEmissions()
    self.CalculateConcentrations()

  @property
  def averageWindSpeed(self):
    return self.__averageWindSpeed

  @averageWindSpeed.setter
  def averageWindSpeed(self, v):
    self.__averageWindSpeed = v
    self.CalculateConcentrations()

  @property
  def backgroundConcentration(self):
    return self.__backgroundConcentration

  @backgroundConcentration.setter
  def backgroundConcentration(self, v):
    if v == 'AllZero':
      v = {}
      for Pol in self.emissionFactors.AvailablePollutants:
        v[Pol] = 0
    self.__backgroundConcentration = v
    self.CalculateConcentrations()

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
        self.CalculateConcentrations()

  def GetBounds(self):
    if self.__beingCreated:
      return
    # Find the geographic bounds.
    for index, row in self.geoDataFrame.iterrows():
      if index == 0:
        Bounds_ = reshape(array(row.geometry.bounds), [4, 1])
      else:
        Bounds_ = append(Bounds_, reshape(array(row.geometry.bounds), [4, 1]), 1)
    Bounds = [min(Bounds_[0, :]), min(Bounds_[1, :]), max(Bounds_[2, :]), max(Bounds_[3, :])]
    return Bounds

  def CalculateEmissions(self):
    if self.__beingCreated:
      return
    AvSpeedClasses = self.emissionFactors.AvailableSpeedClasses
    AvPols = self.emissionFactors.AvailablePollutants
    emsAll = {}
    # Remove for loop.


    for index, row in self.geoDataFrame.iterrows():
      if index == 0:

        if row['speedClass'] in AvSpeedClasses:
          speedType = lambda row: row['speedClass']
        elif AvSpeedClasses[0][:2] == 'S_':
          speedType = lambda row: 'S_{:03d}'.format(row['speed'])
        else:
          raise ValueError("The available speed classes in emissionFactors do not agree with the speed class for each road.")
        for Pol in AvPols:
          emsAll[Pol] = []
      ems = self.CalculateEmissionSingle(speedType(row), row['vehicleCounts'], self.emissionFactors)
      for Pol in AvPols:
        emsAll[Pol].append(ems[Pol])


    for Pol in AvPols:
      self.geoDataFrame['emission_{}'.format(Pol)] = emsAll[Pol]

  @staticmethod
  def CalculateEmissionSingle(speedClass, vehicleCounts, emissionFactors):
    Emition = {}
    for Veh in vehicleCounts.vehicleBreakdown:
      vCount = vehicleCounts.counts[Veh]
      for Pol in emissionFactors.AvailablePollutants:

        Emition[Pol] = emissionFactors.Factors[Pol][Veh][speedClass] * vCount
    return Emition

  def CalculateConcentrations(self, distance='RoadEdge', distanceAdd=0, distanceMultiply=1):
    if self.__beingCreated:
      return

    ImpactDistances = [self.__ImpactDistances[s] for s in list(self.geoDataFrame['roadClass'])]
    self.geoDataFrame['impactDistance'] = ImpactDistances
    DispersionCoefficientA = [self.__DispersionCoefficientsAll[s]['A'] for s in list(self.geoDataFrame['roadClass'])]
    self.geoDataFrame['dispersionCoefficient_A'] = DispersionCoefficientA
    DispersionCoefficientB = [self.__DispersionCoefficientsAll[s]['B'] for s in list(self.geoDataFrame['roadClass'])]
    self.geoDataFrame['dispersionCoefficient_B'] = DispersionCoefficientB
    DispersionCoefficientC = [self.__DispersionCoefficientsAll[s]['C'] for s in list(self.geoDataFrame['roadClass'])]
    self.geoDataFrame['dispersionCoefficient_C'] = DispersionCoefficientC
    DispersionCoefficientAlpha = [self.__DispersionCoefficientsAll[s]['Alpha'] for s in list(self.geoDataFrame['roadClass'])]
    self.geoDataFrame['dispersionCoefficient_Alpha'] = DispersionCoefficientAlpha

    DD = self.geoDataFrame['width']/2
    if distance == 'RoadEdge':
      distance = DD
    elif distance == 'RoadCentre':
      distance = zeros_like(DD)
    elif len(distance):
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
    DispFactor = DispFactor * 0.62 * 5 * self.geoDataFrame['treeFactor'] / self.averageWindSpeed

    AvPols = self.emissionFactors.AvailablePollutants
    for Pol in AvPols:
      self.geoDataFrame['concentration_{}'.format(Pol)] = DispFactor*self.geoDataFrame['emission_{}'.format(Pol)] + self.backgroundConcentration[Pol]

  @staticmethod
  def Shape2RoadNetwork(filename):

    # Requires self for input checking.

    # The default road segment is worst case.
    DefaultD = {'roadID': 1, 'roadName': 'Unnamed Road',
                'roadClass': 'NarrowCanyon', 'speedClass': 'Normal',
                'speed': 10, 'width': 5, 'stagnation': 50, 'treeFactor': 1.5,}

    # Define fieldnames which will be looked for within the shape file to
    # assign attributes.
    AttrNames = {}
    AttrNames['roadID'] = ['roadID', 'RoadID', 'OBJECTID']
    AttrNames['roadName'] = ['roadName', 'RoadName', 'ROADNAME']
    AttrNames['roadClass'] = ['roadName', 'Road_Type', 'RoadType', 'Road_Class', 'RoadClass']
    AttrNames['speedClass'] = ['speedClass', 'Speed_Clas', 'SpeedClass', 'Speed_Type', 'SpeedType']
    AttrNames['speed'] = ['speed', 'Speed', 'SPEED']
    AttrNames['width'] = ['WIDTH', 'RoadWidth', 'Width', 'width']
    AttrNames['stagnation'] = ['stagnation', 'Stagnation']
    AttrNames['treeFactor'] = ['treeFactor', 'TreeFactor', 'Tree_Facto']

    # Define fieldnames which will be looked for within the shape file to assign
    # vehicle counts. It may be neccesary to make other breakdowns allowed.
    VehNames = {}
    VehNames['MCycle'] = ['MCycle']
    VehNames['Car'] = ['Car']
    VehNames['LGV'] = ['LGV']
    VehNames['RHGV_2X'] = ['RHGV_2X']
    VehNames['RHGV_3X'] = ['RHGV_3X']
    VehNames['RHGV_4X'] = ['RHGV_4X']
    VehNames['AHGV_34X'] = ['AHGV_34X']
    VehNames['AHGV_5X'] = ['AHGV_5X']
    VehNames['AHGV_6X'] = ['AHGV_6X']
    VehNames['Bus'] = ['Bus']

    # The default road segment is worst case.
    DefaultD = {'roadID': 1, 'roadName': 'Unnamed Road',
                'roadClass': 'NarrowCanyon', 'speedClass': 'Normal',
                'speed': 10, 'width': 5, 'stagnation': 50, 'treeFactor': 1.5}

    # The default vehicle counts object will be empty
    DefaultVC = {'MCycle': 0, 'Car': 0, 'LGV': 0, 'Bus': 0,
                 'RHGV_2X': 0, 'RHGV_3X': 0, 'RHGV_4X': 0,
                 'AHGV_34X': 0, 'AHGV_5X': 0, 'AHGV_6X': 0,}

    # Read the shape file.
    shpFile = gpd.read_file(filename)
    shpKeys = shpFile.keys()

    UseAttrs = {}
    for key, values in AttrNames.items():
      Got = False
      for v in values:
        if v in shpKeys:
          UseAttrs[key] = v
          Got = True
      if not Got:
        print 'No field found for attribute "{}"; default value "{}" will be used.'.format(key, DefaultD[key])

    UseVehs = {}
    for key, values in VehNames.items():
      Got = False
      for v in values:
        if v in shpKeys:
          UseVehs[key] = v
          Got = True
      if not Got:
        print 'No field found for vehicle "{}"; default value of 0 will be used.'.format(key)

    Ds = []
    for index, rd in shpFile.iterrows():
      # For each road in the shape file.
      # Set up a VehicleCounts object.
      VC = DefaultVC
      for VehA, VehB in UseVehs.items():
        VC[VehA] = rd[VehB]
      VC = VehicleCounts(VC)

      # Create a dictionary for this road segment.
      D = DefaultD.copy()
      D['geometry'] = rd.geometry
      D['vehicleCounts'] = VC

      AttKeys = UseAttrs.keys()
      # Error checking.
      for AttrA, AttrB in UseAttrs.items():
        V = rd[AttrB]
        if AttrA == 'roadClass':
          V = V.replace(" ", "")
        D[AttrA] = V
      if 'roadID' in AttKeys:
        rID = rd[UseAttrs['roadID']]
      else:
        rID = index
      D['roadID'] = rID
      # Add the dictionary to a list of dictionaries.
      Ds.append(D)
    # Create a dataframe from the list of dictionaries.
    RN = RoadNetwork(Ds)
    return RN