# -*- coding: utf-8 -*-
"""
Vehicle Count Class
Created on Tue Oct 04 10:22:29 2016

@author: edward.barratt
"""
from numpy import isscalar


class VehicleCounts(object):
  """
  VehicleCounts

  Attributes:
    counts           - dict
    countsScaled     - dict
    duration         - scalar (hours)
    scaling          - dict
    scaling24        - scalar
    vehicleBreakdown - list
  """

  #__VehsToEmissionFactors = {}

  def __init__(self, counts, duration=24, scaling24=1, scaling='Default'):
    self.__beingCreated = True
    self.__vehicleBreakdown = []
    self.__countsScaled = {}
    self.counts = counts
    self.duration = duration
    self.scaling = scaling
    self.scaling24 = scaling24
    self.__beingCreated = False

  @property
  def counts(self):
    return self.__counts

  @counts.setter
  def counts(self, V):
    self.__counts = V
    self.__vehicleBreakdown = V.keys()
    if not self.__beingCreated:
      for veh in self.vehicleBreakdown:
        self.countsScaled[veh] = self.scaling24 * self.scaling * self.counts[veh]

  @property
  def vehicleBreakdown(self):
    return self.__vehicleBreakdown

  @property
  def scaling(self):
    return self.__scaling

  @scaling.setter
  def scaling(self, V):
    scaling = {}
    for veh in self.vehicleBreakdown:
      if V == 'Default':
        scaling[veh] = 1
      elif isinstance(V, dict):
        scaling[veh] = V[veh]
      elif isscalar(V):
        scaling[veh] = V
      else:
        raise ValueError("Scaling value not understood.")
      if not self.__beingCreated:
        self.countsScaled[veh] = scaling[veh] * self.scaling24 * self.counts[veh]
    self.__scaling = scaling

  @property
  def scaling24(self):
    return self.__scaling24

  @scaling24.setter
  def scaling24(self, V):
    self.__scaling24 = V
    for veh in self.vehicleBreakdown:
      self.countsScaled[veh] = V * self.scaling[veh] * self.counts[veh]

  @property
  def countsScaled(self):
    return self.__countsScaled