"""
TestElasticityHighOrder2DPython
is based on the scene 
/cygdrive/c/hdeling/SOFA/SOFA-HighOrder/examples/TetrahedronAssemblyBenchmark/TestElasticityHighOrder2D.xml
but it uses the SofaPython plugin. 
Further informations on the usage of the plugin can be found in 
sofa/applications/plugins/SofaPython/doc/SofaPython.pdf
To lance the scene, type 
runSofa /cygdrive/c/hdeling/SOFA/SOFA-HighOrder/examples/TetrahedronAssemblyBenchmark/TestElasticityHighOrder2DPython.scn

The current file has been written by the python script
../../../fork-github/src/sofa/applications/plugins/SofaPython/scn2python.py
Author of scn2python.py: Christoph PAULUS, christoph.paulus@inria.fr
"""

import Sofa

class TestElasticityHighOrder2D (Sofa.PythonScriptController):

    def createGraph(self,rootNode):

        # rootNode
        rootNode.createObject('RequiredPlugin', pluginName='SofaHighOrderTopology', name='SofaHighOrderTopology plugin')
        rootNode.createObject('RequiredPlugin', pluginName='SofaHighOrderFEM', name='SofaHighOrderFEM plugin')
        rootNode.createObject('RequiredPlugin', pluginName='Electrophysiology')
        rootNode.createObject('Object', shadows='0', type='LightManager', ambient='1 1 1 1')

        # rootNode/MecaQuadNode
        MecaQuadNode = rootNode.createChild('MecaQuadNode')
        MecaQuadNode.gravity = '0 0 0'
        MecaQuadNode.createObject('MechanicalObject', name='mecaObj', template='Vec2d', tags='meca')
        MecaQuadNode.createObject('CubeTopology', nx='10', ny='10', nz='1', name='grid', min='0 0 0', max='2 2 0', tags='meca')
        MecaQuadNode.createObject('QuadSetTopologyContainer', quads='@grid.quads', edges='@grid.edges', name='quadContainer')
        MecaQuadNode.createObject('QuadSetGeometryAlgorithms', template='Vec2d')

        # rootNode/MecaQuadNode/MecaTriangleNode
        MecaTriangleNode = MecaQuadNode.createChild('MecaTriangleNode')
        MecaTriangleNode.createObject('TriangleSetTopologyContainer', position='@../mecaObj.position', name='triangleContainer')
        MecaTriangleNode.createObject('TriangleSetTopologyModifier')
        MecaTriangleNode.createObject('Quad2TriangleTopologicalMapping', input='@../quadContainer', output='@triangleContainer')
        MecaTriangleNode.createObject('TriangleSetGeometryAlgorithms', drawEdges='1', template='Vec2d')

        # rootNode/MecaQuadNode/MecaTriangleNode/HighOrder
        HighOrder = MecaTriangleNode.createChild('HighOrder')
        HighOrder.createObject('CGLinearSolver', threshold='1.0e-9', tolerance='1.0e-9', name='linear solver', iterations='500')
        HighOrder.createObject('StaticSolver', applyIncrementFactor='1')
        HighOrder.createObject('TriangleSetTopologyModifier')
        HighOrder.createObject('HighOrderTriangleSetTopologyContainer', name='ContainerBezier')
        HighOrder.createObject('Mesh2HighOrderTopologicalMapping', bezierTriangleDegree='2', output='@ContainerBezier', name='topoMapping', input='@../triangleContainer')
        HighOrder.createObject('MechanicalObject', showObject='0', name='BezierMeca', template='Vec2d')
        HighOrder.createObject('BezierTriangleSetGeometryAlgorithms', drawEdges='true', drawControlPoints='true', name='GeomAlgo', template='Vec2d', drawColorEdges='1 0 0')
        HighOrder.createObject('HighOrderTriangularCorotationalFEMForceField', integrationMethod='standard', name='Elasticity', integrationOrder='4', youngModulus='4', printLog='1', template='Vec2d', numericalIntegrationMethod='3', method='linear', poissonRatio='0.45')
        HighOrder.createObject('BoxROI', box='-0.01 -0.01 -0.01 2 0.01 0.01', drawBoxes='true', name='ROIZL', template='Vec3d')
        HighOrder.createObject('FixedConstraint', indices='@ROIZL.indices')

        return 0;

    def onKeyPressed(self, c):
        return 0;

    def onKeyReleased(self, c):
        return 0;

    def onLoaded(self, node):
        return 0;

    def onMouseButtonLeft(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseButtonRight(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseButtonMiddle(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseWheel(self, mouseX,mouseY,wheelDelta):
        return 0;

    def onGUIEvent(self, strControlID,valueName,strValue):
        return 0;

    def onBeginAnimationStep(self, deltaTime):
        return 0;

    def onEndAnimationStep(self, deltaTime):
        return 0;

    def onScriptEvent(self, senderNode, eventName,data):
        return 0;

    def initGraph(self, node):
        return 0;

    def bwdInitGraph(self, node):
        return 0;

    def storeResetState(self):
        return 0;

    def reset(self):
        return 0;

    def cleanup(self):
        return 0;
