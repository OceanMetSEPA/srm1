# -*- coding: utf-8 -*-
"""
EmissionFactorClass
Created on Mon Oct 03 14:35:18 2016

@author: edward.barratt
"""
import os

class EmissionFactorClass(object):
  """
  EmissionFactorClass

  Attributes:
    AvailablePollutants
    AvailableVehicles
    AvailableSpeedClasses
    Name
    Units
    Factors   - dict
  """

  def __init__(self, Factors='NoneSet', Name='Unnamed', Units='g/km/vehicle'):


    # Pre define some hidden variables which will become populated when factors
    # are imported.
    self.__AvailablePollutants = []
    self.__AvailableVehicles = []
    self.__AvailableSpeedClasses = []

    # Get the factors
    if Factors == 'NoneSet':
      if Name == 'Unnamed':
        raise ValueError("Either Factors should be defined, as a dictionary, ",
                         "or an existing emission factor class name should be defined.")
      else:
        [Blah, self.Factors, Units] = PreSetEmissionFactorClasses(Name)
        self.__Name = Name
    else:
      self.Factors = Factors;
      self.Name = Name;
    self.Units = Units

  @property
  def Name(self):
    return self.__Name

  @Name.setter
  def Name(self, V):
    try:
      [FName, Blah] = PreSetEmissionFactorClasses(V)
    except ValueError:
      self.__Name = V
    else:
      print 'a'
      print ("A pre-defined emission factor class with that name is already "
             "saved at {}, are you sure you want to create another emission "
             "factor class with that name?").format(FName)
      # Add question prompt

  @property
  def Factors(self):
    return self.__Factors

  @Factors.setter
  def Factors(self, V):
    self.__AvailablePollutants = V.keys()
    self.__AvailableVehicles = V[self.AvailablePollutants[0]].keys()
    self.__AvailableSpeedClasses = V[self.AvailablePollutants[0]][self.AvailableVehicles[0]].keys()
    self.__Factors = V

  @property
  def AvailablePollutants(self):
    return self.__AvailablePollutants

  @property
  def AvailableVehicles(self):
    return self.__AvailableVehicles

  @property
  def AvailableSpeedClasses(self):
    return self.__AvailableSpeedClasses

  # End of class EmissionFactorClass:


def PreSetEmissionFactorClasses(Name):
  FileName = 'Data\EF_{}.txt'.format(Name)
  Units='g/km/vehicle'
  Factors = {}
  if not os.path.isfile(FileName):
    raise ValueError("No pre set emission factor file can be found for assigned name.")
  with open(FileName, 'r') as f:
    lines = f.readlines()
    Data = False
    for l in lines:
      l = l.strip()
      if l == 'DataStart':
        Data = True
      elif l == 'DataEnd':
        Data = False
      elif l[:6] == 'Units:':
        Units = l[6:].strip()
      elif Data:
        if l[:14] == 'PollutantName:':
          PName = l[14:].strip()
          Factors[PName] = {}
        elif l[:17] == 'VehicleClassName:':
          VName = l[17:].strip()
          Factors[PName][VName] = {}
        else:
          SC = l.split(': ')
          Factors[PName][VName][SC[0]] = float(SC[1])
  return FileName, Factors, Units

