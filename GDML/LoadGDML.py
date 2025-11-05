import pyg4ometry

r = pyg4ometry.gdml.Reader("/l/triton_work/RadEx/RadFETs/RadFETs.gdml")
l = r.getRegistry().getWorldVolume()
v = pyg4ometry.visualisation.VtkViewerColouredMaterial()
v.addLogicalVolume(l)
v.view()