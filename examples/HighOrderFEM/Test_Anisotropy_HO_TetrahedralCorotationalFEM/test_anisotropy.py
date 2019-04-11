# -*- coding: utf-8 -*-
import os
import sys

#   STLIB IMPORT
from stlib.scene import Scene
from stlib.solver import DefaultSolver
from stlib.physics.deformable import ElasticMaterialObject
from stlib.physics.constraints import FixedBox

mesh = 'myBeam.vtk'

def highOrderCube(node,poissonRatio,youngModulus,totalMass,name="highOrderCube",translation=[0,0,0],rotation=[0,0,0],
    numericalIntegrationMethod='Tetrahedron Gauss',
    forceAffineAssemblyForAffineElements= False, # if true affine tetrahedra are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod
    integrationMethod='analytical', # analytical # numerical # standard !! if degree ==  1 --> force to AFFINE_ELEMENT_INTEGRATION
    method="qr", # linear # qr # polar # polar2 # none
    integrationOrder=1,
    oneRotationPerIntegrationPoint= False,
    withAnysotropy=True,
    elasticitySymmetry='transverseIsotropic', # \"isotropic\"  or \"transverseIsotropic\" or \"orthotropic\" or \"cubic\
    youngModulusLongitudinal = None,
    anisotropyParameters=[1000,0.45,300], # [youngModulusLongitudinal, poissonRatioTransverseLongitudinal, shearModulusTransverse]
    anisotropyDirections = [[1,0,0]] # vector<Coord>
    ):

    print("highOrderCube : "+name)
    print("    - youngModulus : " + str(youngModulus))
    print("    - youngModulusLongitudinal : " + str(youngModulusLongitudinal))
    print("    - poissonRatio : " + str(poissonRatio))
    
    if youngModulusLongitudinal != None: 
        # http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
        # 3.2.14 Stress-strain relations for linear elastic Transversely Isotropic Material
        poissonRatioTransverseLongitudinal = (poissonRatio*youngModulusLongitudinal)/youngModulus
        shearModulusTransverse = youngModulusLongitudinal/(2*(1+poissonRatioTransverseLongitudinal))

        anisotropyParameters = [youngModulusLongitudinal,poissonRatioTransverseLongitudinal,shearModulusTransverse]

        print("    - poissonRatioTransverseLongitudinal : " + str(poissonRatioTransverseLongitudinal))
        print("    - shearModulusTransverse : " + str(shearModulusTransverse))

    # MESH TOPO
    highOrderCube = node.createChild(name)
    gc = highOrderCube.createObject('MeshVTKLoader',name="mesh",filename=mesh,translation=translation,rotation=rotation)
    tetstc = highOrderCube.createObject('TetrahedronSetTopologyContainer',name="Container1",tetrahedra="@mesh.tetrahedra",position="@mesh.position")
    tetsga = highOrderCube.createObject('TetrahedronSetGeometryAlgorithms',name="GeomAlgo")    
    mecaTet = highOrderCube.createObject('MechanicalObject',name="dofs",rotation=rotation,translation=translation)

    highOrder_node = highOrderCube.createChild('highOrder_node')
    DefaultSolver(highOrder_node)


    highOrder_node.createObject('HighOrderTetrahedronSetTopologyContainer',name="ContainerBezier")
    highOrder_node.createObject('Mesh2HighOrderTopologicalMapping',name="topoMapping",input="@"+name+"/Container1",output="@ContainerBezier",
                                            bezierTetrahedronDegree=integrationOrder)

    highOrder_node.createObject('MechanicalObject',name="hodofs")
    highOrder_node.createObject('UniformMass', totalMass=totalMass)

    highOrder_node.createObject('LagrangeTetrahedronSetGeometryAlgorithms',name="GeomAlgo",
                                            drawControlPoints=1,
                                            drawEdges=1,
                                            drawColorEdges="1 0 0")

    high_fem = highOrder_node.createObject('HighOrderTetrahedralCorotationalFEMForceField',
                                            name="Elasticity",  printLog=1,  poissonRatio=poissonRatio,  youngModulus=youngModulus,
                                            integrationMethod=integrationMethod,
                                            method=method,
                                            forceAffineAssemblyForAffineElements=forceAffineAssemblyForAffineElements,
                                            oneRotationPerIntegrationPoint= oneRotationPerIntegrationPoint,
                                            numericalIntegrationMethod=numericalIntegrationMethod,
                                            integrationOrder=integrationOrder)

    if withAnysotropy :
        high_fem.elasticitySymmetry   = elasticitySymmetry
        high_fem.anisotropyParameters = anisotropyParameters
        high_fem.anisotropyDirections = anisotropyDirections
    
    fixedBoxPos=[translation[0], translation[1], translation[2],translation[0]+2, translation[1]+ 15, translation[2]+15]
    FixedBox(   atPositions=fixedBoxPos,
                applyTo=highOrder_node,
                doVisualization=True) 


def normalCube(node,mesh,poissonRatio,youngModulus,totalMass,translation=[40,0,0],rotation=[0,0,0]):
    
    normal_cube = node.createChild('normal')
    DefaultSolver(normal_cube)

    obj= ElasticMaterialObject(
            attachedTo=normal_cube,
            volumeMeshFileName=mesh,
            name='modelNode',
            translation=translation,
            rotation=rotation,
            totalMass=totalMass,
            withConstrain=False,
            # surfaceMeshFileName='myCube.obj', 
            surfaceColor=[0.0, 0.7, 0.7, 0.5],
            poissonRatio=poissonRatio,
            youngModulus=youngModulus)

    fixedBoxPos=[translation[0], translation[1], translation[2],translation[0]+2, translation[1]+ 15, translation[2]+15]
    FixedBox(   atPositions=fixedBoxPos,
                applyTo=obj,
                doVisualization=True) 


def createScene(rootNode):

    scene = Scene(rootNode, gravity=[0.0,0.0,-9810],
                            dt=0.01,
                            plugins=['SofaHighOrderTopology','SofaHighOrderFEM'])

    scene.VisualStyle.displayFlags="showBehaviorModels showForceFields"
  
    poissonRatio=0.45
    youngModulus=200
    totalMass=0.02
    
    # NORMAL CUBE
    normalCube(rootNode,mesh,poissonRatio,youngModulus,totalMass,translation=[0,0,0])

    # HIGH ORDER Without ANYSOTROPY
    highOrderCube(rootNode,poissonRatio,youngModulus,totalMass,name="HO_isotropic",elasticitySymmetry='isotropic',
                translation=[0,30,0],youngModulusLongitudinal=200,withAnysotropy=False) #,integrationMethod='standard')

    # HIGH ORDER With ANYSOTROPY CUBIC
    highOrderCube(rootNode,poissonRatio,youngModulus,totalMass,name="HO_with_ani_cubic",elasticitySymmetry='cubic',
                translation=[0,60,0],anisotropyParameters=[1])

    # HIGH ORDER With ANYSOTROPY TRANSVERSE
    highOrderCube(rootNode,poissonRatio,youngModulus,totalMass,name="HO_with_ani_transverse",elasticitySymmetry='transverseIsotropic',
                translation=[0,90,0],youngModulusLongitudinal=200)

    # # HIGH ORDER 
    # highOrderCube(rootNode,poissonRatio,youngModulus,totalMass,name="highOrder_with_anisotropy_diff",
    #             translation=[0,120,0],youngModulusLongitudinal=100)

    # # HIGH ORDER 
    # highOrderCube(rootNode,poissonRatio,youngModulus,totalMass,name="highOrderCube_4",
    #             translation=[0,120,0],youngModulusLongitudinal=100)

