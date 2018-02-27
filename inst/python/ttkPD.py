#### import the simple module from the paraview
from paraview.simple import *

if __name__ == '__main__':
    ### Inputs
    inputs = sys.argv[1:]

    #path = '/home/duynguyen/ResearchLocal/HiCDA/Simulations/DeriveData/simulation/'
    path = inputs[0]
    # create a new 'Legacy VTK Reader'
    #simvtk = LegacyVTKReader(FileNames=['/home/duynguyen/Dropbox/Research/HiCDA/Simulations/DerivedData/simulation/sim.vtk'])
    simvtk = LegacyVTKReader(FileNames=[path + 'sim.vtk'])


    # create a new 'TTK PersistenceDiagram'
    tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=simvtk)

    # Properties modified on tTKPersistenceDiagram1
    tTKPersistenceDiagram1.InputOffsetField = ''

    # TODO save data
    SaveData(path + 'PD-CellData.csv', proxy=tTKPersistenceDiagram1, FieldAssociation='Cells')

    ## create a new 'Threshold'
    #threshold1 = Threshold(Input=tTKPersistenceDiagram1)
    #
    ## Properties modified on threshold1
    #threshold1.Scalars = ['CELLS', 'Persistence']
    #threshold1.ThresholdRange = [0.0, 1.0]
    #
    #
    ## create a new 'TTK Merge and Contour Tree (FTM)'
    #tTKMergeandContourTreeFTM1 = TTKMergeandContourTreeFTM(Input=threshold1)
    #
    ## get active source.
    #tTKMergeandContourTreeFTM1_1 = GetActiveSource()
    #
    ## save data
    #SaveData(path+'thresholdGlobal.csv', proxy=OutputPort(tTKMergeandContourTreeFTM1_1, 2))
    #
