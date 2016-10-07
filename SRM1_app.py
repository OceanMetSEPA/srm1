# -*- coding: utf-8 -*-
"""
Created on Thu Oct 06 16:02:32 2016

@author: edward.barratt
"""

from flask import Flask, render_template, redirect, request
from SRM1RoadNetwork import RoadNetwork
import folium
from numpy import mean

app = Flask(__name__)
app.ShapeFile = 'Data\151119AbdnRoadsForDutchModel.shp'
app.RoadNetwork = RoadNetwork.Shape2RoadNetwork(ShapeFileName)


def PlotNetwork():

  LatMin = app.RoadNetwork.Bounds[0]
  LonMin = app.RoadNetwork.Bounds[1]
  LatMax = app.RoadNetwork.Bounds[2]
  LonMax = app.RoadNetwork.Bounds[3]
  LatMean = mean([LatMin, LatMax])
  LonMean = mean([LonMin, LonMax])

  # Create a map
  map_osm = folium.Map(zoom_start = 6, location=[LatMean, LonMean])

  # Create a marker cluster
  marker_cluster = folium.MarkerCluster("Sites").add_to(map_osm)
  # Add a marker for each site.
  for siteI in xrange(len(SiteInfo.index)):
    if siteI >= 400:
      break
    Lat = SiteInfo['latitude'][siteI]
    Lon = SiteInfo['longitude'][siteI]
    folium.Marker([Lat, Lon]).add_to(marker_cluster)
  map_osm.fit_bounds([[LatMin, LonMin], [LatMax, LonMax]])

  #TempHTMLFile = 'templates/TempMap_{}.html'.format(RandomString())
  TempHTMLFile = 'templates/Map.html'
  map_osm.save(TempHTMLFile)
  app.MapFile = TempHTMLFile

    # Tool tip
    #hover = HoverTool(tooltips=[("Site ID", '@id'),
    #                            ("name", '@name'),
    #                            ("Longitude", '@longitude'),
    #                            ("Latitude", '@latitude'),
    #                            ("County", '@unitaryAuthArea'),
    #                            ("region", '@region')])
    p = figure(width=700, height=700, title="SRM1 Road Network",
               tools=['pan', 'box_zoom', 'wheel_zoom', 'save', 'reset', hover])
    p.circle('longitude', 'latitude', size=7, color="firebrick", alpha=0.5, source=source)
    script, div = components(p)
    return script, div



@app.route('/')
def main():
  return redirect('/index')

@app.route('/index')
def index():
  #RegionDropDown = CreateRegionDropDown()
  SitePlotScript, SitePlotDiv = PlotSites()
  return render_template('index.html', CountyForm=RegionDropDown,
                                       NetworkPlotScript=SitePlotScript,
                                       NetworkPlotDiv=SitePlotDiv)


if __name__ == '__main__':
  app.run(debug=True)
