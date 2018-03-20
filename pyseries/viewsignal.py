"""
Script for signal visualization
"""
import RDCsignals as rdc
import RDClectura as lec
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
#
def mouseMoved(evt):
    pos = evt[0]
    mousePoint = vb.mapSceneToView(pos)
    vLine.setPos(mousePoint.x())
    hLine.setPos(mousePoint.y())
    label.setText("<span style='font-size: 14pt; color: red'> x = %0.2f, <span style='color: green'> y = %0.2f</span>" % (mousePoint.x(), mousePoint.y()))
#
def update():
    region.setZValue(10)
    minX, maxX = region.getRegion()
    p1.setXRange(minX, maxX, padding=0)
#
def updateRegion(window, viewRange):
    rgn = viewRange[0]
    region.setRegion(rgn)
"""
Read the time series
"""
ndatos , aa = lec.readsignal_VEC('pulso' , 1.0)
#ndatos , aa = lec.readsignal_OPENSEES('accele' , 1.0)
n = len(aa)
y = rdc.linebaseace(aa[:])
x=np.zeros([n], dtype=float)
x=np.arange(0,n*0.01,0.01)
"""
Generate layout
"""
win = pg.GraphicsWindow()
win.setWindowTitle('RDC_series')
label = pg.LabelItem(justify = "right")
win.addItem(label)
p1 = win.addPlot(row=1, col=0)
p2 = win.addPlot(row=2, col=0)
p1.plot(x , y , pen="r")
p2.plot(x , y , pen="w")
"""
Add movable region
"""
region = pg.LinearRegionItem()
region.setZValue(10)
p2.addItem(region, ignoreBounds=True)
p1.setAutoVisible(y=True)
#
region.sigRegionChanged.connect(update)
p1.sigRangeChanged.connect(updateRegion)
region.setRegion([20, 40])
"""
Generate the crosshair
"""
vLine = pg.InfiniteLine(angle=90, movable=True)
hLine = pg.InfiniteLine(angle= 0, movable=True)
p1.addItem(vLine, ignoreBounds=True)
p1.addItem(hLine, ignoreBounds=True)
p1.showGrid(1 , 1 , 1.0)
p2.showGrid(1 , 1 , 0.5)
#
vb = p1.vb
proxy = pg.SignalProxy(p1.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)
QtGui.QApplication.instance().exec_()
pg.QtGui.QApplication.processEvents()