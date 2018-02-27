import sys
import io
import getopt
from paraview.simple import *
#### import the simple module from the paraview
from paraview.simple import *



if __name__ == '__main__':
    ### Inputs
    inputs = sys.argv[1:]
    # path = '/home/duynguyen/ResearchLocal/HiCDA/Simulations/DeriveData/simulation/'
    path = inputs[0]
    # create a new 'Legacy VTK Reader'
    simvtk = LegacyVTKReader(FileNames=[path + 'sim.vtk'])

    # create a new 'TTK PersistenceDiagram'
    tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=simvtk)

    # Properties modified on tTKPersistenceDiagram1
    tTKPersistenceDiagram1.InputOffsetField = ''

    # create a new 'Threshold'
    threshold1 = Threshold(Input=tTKPersistenceDiagram1)

    # Properties modified on threshold1
    threshold1.Scalars = ['CELLS', 'Persistence']
    threshold_range = [float(inputs[1]), float(inputs[2])]
    threshold1.ThresholdRange = threshold_range

    # create a new 'TTK TopologicalSimplification'
    tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=simvtk,
        Constraints=threshold1)

    # Properties modified on tTKTopologicalSimplification1
    tTKTopologicalSimplification1.InputOffsetField = ''

    # create a new 'TTK MorseSmaleComplex'
    tTKMorseSmaleComplex1 = TTKMorseSmaleComplex(Input=tTKTopologicalSimplification1)

    # get active source.
    tTKMorseSmaleComplex1_1 = GetActiveSource()

    # save data
    #SaveData('/home/duynguyen/ResearchLocal/HiCDA/Simulations/DeriveData/simulation/_MSComplex.csv', proxy=OutputPort(tTKMorseSmaleComplex1_1, 3))
    # save data
    #fileout = 'MSComplex' + str(threshold_range[0]) + '-' + str(threshold_range[1])
    fileout = 'MSComplexPLevel1'
    SaveData(path+fileout+'.csv', proxy=OutputPort(tTKMorseSmaleComplex1_1, 3))
