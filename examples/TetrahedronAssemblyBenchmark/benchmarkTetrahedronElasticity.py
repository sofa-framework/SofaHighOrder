import Sofa
import sys

############################################################################################
# this is a PythonScriptController example script
############################################################################################




############################################################################################
# following defs are used later in the script
############################################################################################

ffList= []
topoListTetra= []
topoListTrian= []
filename= "dummy"

def highOrderElasticityTrianScene(node,counter,degree,resolution,isForceField,integrationType,numericalIntegrationOrder,lumped=0):
#	print "highOrderElasticityTrianScene"
	meshTopology = node.createChild("MeshTopology" +str(counter))
	res=2*resolution
	resStr=str(res)+" "+str(res)+" 0"

	cubeTopo = meshTopology.createObject('GenerateGrid',name="grid",template="Vec2d",resolution=resStr,min="0 0 0",max="2 2 0",tags="meca")
	trianstc = meshTopology.createObject('TriangleSetTopologyContainer',name="triangleContainer",position="@grid.output_position",triangles="@grid.triangles")
	mecaTrian = meshTopology.createObject('MechanicalObject',name="mecaObj",template='Vec2d',tags="meca")
	triansga = meshTopology.createObject('TriangleSetGeometryAlgorithms',template="Vec2d")	
	
	# rootNode/MecaQuadNode/MecaTriangleNode/HighOrder
	HighOrder = meshTopology.createChild('HighOrder')
	HighOrder.createObject('TriangleSetTopologyModifier')
	hotetstc=HighOrder.createObject('HighOrderTriangleSetTopologyContainer', name='ContainerBezier')
	topoListTrian.append(hotetstc)
	HighOrder.createObject('Mesh2HighOrderTopologicalMapping', bezierTriangleDegree=str(degree), output='@ContainerBezier', name='topoMapping', input='@../triangleContainer')
	HighOrder.createObject('MechanicalObject', showObject='0', name='BezierMeca', template='Vec2d')
	HighOrder.createObject('LagrangeTriangleSetGeometryAlgorithms', drawEdges='true', drawControlPoints='true', name='GeomAlgo', template='Vec2d', drawColorEdges='1 0 0')
	if (isForceField):
		hotcfemff=HighOrder.createObject('HighOrderTriangularCorotationalFEMForceField', name='Elasticity', integrationMethod=str(integrationType),forceAffineAssemblyForAffineElements="0",numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder), youngModulus='4', printLog='1', template='Vec2d', method='linear', poissonRatio='0.45')
	else:
		hotcfemff = HighOrder.createObject('HighOrderMeshMatrixMass',name="Mass",template='Vec2d',printLog="1",lumping=str(lumped),massDensity="1.0",forceAffineAssemblyForAffineElements="0",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder))

	ffList.append(hotcfemff)  
	return 0
def highOrderDiffusionTrianScene(node,counter,degree,resolution,anisotropic,integrationType,numericalIntegrationOrder):
#	print "highOrderDiffusionTrianScene"
	meshTopology = node.createChild("MeshTopology" +str(counter))
	res=2*resolution
	resStr=str(res)+" "+str(res)+" 0"

	cubeTopo = meshTopology.createObject('GenerateGrid',name="grid",template="Vec2d",resolution=resStr,min="0 0 0",max="2 2 0",tags="meca")
	trianstc = meshTopology.createObject('TriangleSetTopologyContainer',name="triangleContainer",position="@grid.output_position",triangles="@grid.triangles")
	mecaTrian = meshTopology.createObject('MechanicalObject',name="mecaObj",template='Vec2d',tags="meca")
	triansga = meshTopology.createObject('TriangleSetGeometryAlgorithms',template="Vec2d")	
	
 

	# rootNode/MecaQuadNode/MecaTriangleNode/HighOrder
	HighOrder = meshTopology.createChild('HighOrder')
	HighOrder.createObject('TriangleSetTopologyModifier')
	hotetstc=HighOrder.createObject('HighOrderTriangleSetTopologyContainer', name='ContainerBezier')
	topoListTrian.append(hotetstc)
	HighOrder.createObject('Mesh2HighOrderTopologicalMapping', bezierTriangleDegree=str(degree), output='@ContainerBezier', name='topoMapping', input='@../triangleContainer')
	HighOrder.createObject('MechanicalObject', showObject='0', name='BezierMeca', template='Vec2d')
	HighOrder.createObject('LagrangeTriangleSetGeometryAlgorithms', drawEdges='true', drawControlPoints='true', name='GeomAlgo', template='Vec2d', drawColorEdges='1 0 0')
	highOrder1D_node = HighOrder.createChild('elec_node')
	elecHoTet = highOrder1D_node.createObject('MechanicalObject',name="elecDofs",template="Vec1d",position="0",tags="EPModel")
	btsga2 = highOrder1D_node.createObject('LagrangeTriangleSetGeometryAlgorithms',template="Vec1d")
	if (anisotropic):
		hotdfemff = highOrder1D_node.createObject('HighOrderTriangularDiffusionForceField',name="Diffusion",forceAffineAssemblyForAffineElements="0",anisotropy="isotropy",anisotropyDirection="1 0 0",anisotropyParameters="3",template="Vec1d",printLog="1",diffusivity="0.01",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder),tags="EPModel" )
	else:
		hotdfemff = highOrder1D_node.createObject('HighOrderTriangularDiffusionForceField',name="Diffusion",template="Vec1d",printLog="1",forceAffineAssemblyForAffineElements="0",diffusivity="0.01",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder),tags="EPModel" )
	ffList.append(hotdfemff)
	return 0	

def highOrderElasticityTetraScene(node,counter,degree,resolution,isForceField,integrationType,numericalIntegrationOrder,lumped=0):
#	print "highOrderElasticityTetraScene"
	meshTopology = node.createChild("MeshTopology" +str(counter))
	resHeightvalue=1*resolution
	resCircumferentialValue=2*resolution
	resRadialValue=1*resolution
	gc = meshTopology.createObject('GenerateCylinder',name="Cylinder",radius="0.2",height="1.0",resHeight=str(resHeightvalue),resCircumferential=str(resCircumferentialValue),resRadial=str(resRadialValue))
	tetstc = meshTopology.createObject('TetrahedronSetTopologyContainer',name="Container1",tetrahedra="@Cylinder.tetrahedra",position="@Cylinder.output_TetrahedraPosition")
	tetsga = meshTopology.createObject('TetrahedronSetGeometryAlgorithms',name="GeomAlgo",drawEdges="0")	
	mecaTet = meshTopology.createObject('MechanicalObject',name="dofs")
	highOrder_node = meshTopology.createChild('particle_node')
	hotetstc = highOrder_node.createObject('HighOrderTetrahedronSetTopologyContainer',name="ContainerBezier")
	topoListTetra.append(hotetstc)
	m2hotm = highOrder_node.createObject('Mesh2HighOrderTopologicalMapping',name="topoMapping",input="@../Container1",output="@ContainerBezier",bezierTetrahedronDegree=str(degree))
	mecaHoTet = highOrder_node.createObject('MechanicalObject',name="hodofs")
	btsga = highOrder_node.createObject('LagrangeTetrahedronSetGeometryAlgorithms',name="GeomAlgo",drawControlPoints="true",drawEdges="true",drawColorEdges="1 0 0")
	if (isForceField):
		hotcfemff = highOrder_node.createObject('HighOrderTetrahedralCorotationalFEMForceField',name="Elasticity",printLog="1",forceAffineAssemblyForAffineElements="0",poissonRatio="0.45",youngModulus="4",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder),method="linear"  )
	else:
		hotcfemff = highOrder_node.createObject('HighOrderMeshMatrixMass',name="Mass",printLog="1",lumping=str(lumped),massDensity="1.0",forceAffineAssemblyForAffineElements="0",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder))
	ffList.append(hotcfemff)
	return 0

def highOrderDiffusionTetraScene(node,counter,degree,resolution,anisotropic,integrationType,numericalIntegrationOrder):
#	print "highOrderDiffusionTetraScene"
	meshTopology = node.createChild("MeshTopology" +str(counter))
	resHeightvalue=1*resolution
	resCircumferentialValue=2*resolution
	resRadialValue=1*resolution
	gc = meshTopology.createObject('GenerateCylinder',name="Cylinder",radius="0.2",height="1.0",resHeight=str(resHeightvalue),resCircumferential=str(resCircumferentialValue),resRadial=str(resRadialValue))
	tetstc = meshTopology.createObject('TetrahedronSetTopologyContainer',name="Container1",tetrahedra="@Cylinder.tetrahedra",position="@Cylinder.output_TetrahedraPosition")
	tetsga = meshTopology.createObject('TetrahedronSetGeometryAlgorithms',name="GeomAlgo",drawEdges="0")	
	mecaTet = meshTopology.createObject('MechanicalObject',name="dofs")
	highOrder_node = meshTopology.createChild('particle_node')
	hotetstc = highOrder_node.createObject('HighOrderTetrahedronSetTopologyContainer',name="ContainerBezier")
	topoListTetra.append(hotetstc)
	m2hotm = highOrder_node.createObject('Mesh2HighOrderTopologicalMapping',name="topoMapping",input="@../Container1",output="@ContainerBezier",bezierTetrahedronDegree=str(degree))
	mecaHoTet = highOrder_node.createObject('MechanicalObject',name="hodofs")
	btsga = highOrder_node.createObject('LagrangeTetrahedronSetGeometryAlgorithms',name="GeomAlgo",drawControlPoints="true",drawEdges="true",drawColorEdges="1 0 0")
	highOrder1D_node = highOrder_node.createChild('elec_node')
	elecHoTet = highOrder1D_node.createObject('MechanicalObject',name="elecDofs",template="Vec1d",position="0",tags="EPModel")
	btsga2 = highOrder1D_node.createObject('LagrangeTetrahedronSetGeometryAlgorithms',template="Vec1d")
	if (anisotropic):
		hotdfemff = highOrder1D_node.createObject('HighOrderTetrahedralDiffusionForceField',name="Diffusion",anisotropy="isotropy",forceAffineAssemblyForAffineElements="0",anisotropyDirection="1 0 0",anisotropyParameters="3",template="Vec1d",printLog="1",diffusivity="0.01",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder),tags="EPModel" )
	else:
		hotdfemff = highOrder1D_node.createObject('HighOrderTetrahedralDiffusionForceField',name="Diffusion",template="Vec1d",printLog="1",forceAffineAssemblyForAffineElements="0",diffusivity="0.01",integrationMethod=str(integrationType),numericalIntegrationMethod="3",integrationOrder=str(numericalIntegrationOrder),tags="EPModel" )
	ffList.append(hotdfemff)
	return 0

############################################################################################
# following defs are optionnal entry points, called by the PythonScriptController component;
############################################################################################


class BenchmarkTetrahedronElasticity(Sofa.PythonScriptController):
	# called once the script is loaded
	def onLoaded(self,node):
		self.counter = 0
		print 'Controller script loaded from node %s'%node.findData('name').value
		return 0

	# optionnally, script can create a graph...
	def createGraph(self,node):
		print 'createGraph called (python side)'

		# first take input from variables array
		variables = self.findData('variables').value
		if len(variables)==6:
			print 'initialization from input file benchmark.params'
			nbtrees= int(variables[1][0])
			degree=int(variables[0][0])
			resolution=int(variables[2][0])
			componentType=variables[3][0]
			integrationType=variables[4][0]
			numericalIntegrationOrder=int(variables[5][0])
		else:	
			print 'initialization from input file benchmark.params'
			f=open("benchmark.params",'r')
			myline=f.readline()
			input=myline.split()
#		print 'At initialization, input file content  ', input[0]
			f.close()
			nbtrees= int(input[1])
			degree=int(input[0])
			resolution=int(input[2])
			componentType=input[3]
			integrationType=input[4]
			numericalIntegrationOrder=int(input[5])
			

		print 'Number of trees= ', nbtrees
		global filename
		filename="benchmarkTetrahedron-"+str(degree)+"-"+str(nbtrees).zfill(2)+"-"+str(resolution).zfill(2)+"-"+componentType+"-"+integrationType+"-"+str(numericalIntegrationOrder)+".bench"
		
		# highOrderElasticityScene(node,1,degree,resolution,componentType,integrationType,numericalIntegrationOrder)
		for x in range(0,nbtrees):
			if componentType=='elasticityTetra':
				highOrderElasticityTetraScene(node,(x),degree,resolution,True,integrationType,numericalIntegrationOrder)
			elif componentType=='massTetra':
				highOrderElasticityTetraScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder)	
			elif componentType=='massTetraLumped':
				highOrderElasticityTetraScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder,1)					
			elif componentType=='diffusionTetra':
				highOrderDiffusionTetraScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder)
			elif componentType=='diffusionTetraAniso':
				highOrderDiffusionTetraScene(node,(x),degree,resolution,True,integrationType,numericalIntegrationOrder)				
			elif componentType=='elasticityTrian':
				highOrderElasticityTrianScene(node,(x),degree,resolution,True,integrationType,numericalIntegrationOrder)
			elif componentType=='massTrian':
				highOrderElasticityTrianScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder)	
			elif componentType=='massTrianLumped':
				highOrderElasticityTrianScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder,1)					
			elif componentType=='diffusionTrian':
				highOrderDiffusionTrianScene(node,(x),degree,resolution,False,integrationType,numericalIntegrationOrder)	
			elif componentType=='diffusionTrianAniso':
				highOrderDiffusionTrianScene(node,(x),degree,resolution,True,integrationType,numericalIntegrationOrder)				
			else:
				print "Unknown component Type",componentType
		try:
			sys.stdout.flush()
		except IOError:
			pass			
		return 0



	# called once graph is created, to init some stuff...
	def initGraph(self,node):
		print 'initGraph called (python side)'
		return 0

	def bwdInitGraph(self,node):
		print 'bwdInitGraph called (python side)'
		assemblyTimeArray=[]
		for ff in ffList:
		   assemblyTimeArray.append(ff.findData('assemblyTime').value)
		print 'mean assembly time=', sum(assemblyTimeArray)/len(assemblyTimeArray), " min=", min(assemblyTimeArray)
		if len(topoListTrian)>0:
			nbElements=len(topoListTrian[0].findData('triangles').value)
			nbPoints=topoListTrian[0].findData('nbPoints').value
		else:
			nbElements=len(topoListTetra[0].findData('tetrahedra').value)
			nbPoints=topoListTetra[0].findData('nbPoints').value	

		print 'nb elements=',nbElements,  'nbPoints=',nbPoints
		try:
			sys.stdout.flush()
		except IOError:
			pass		
		global filename
		filename="benchmarks/"+filename
		f=open(filename,"w")
		f.write( str(sum(assemblyTimeArray)/len(assemblyTimeArray))+" "+ str(min(assemblyTimeArray))+" "+str(nbElements)+" "+str(nbPoints))
		f.close()

		return 0



	def reset(self):
		print 'reset called (python side)'

		return 0

	def cleanup(self):
		print 'cleanup called (python side)'

		return 0


	# called when a GUIEvent is received
	def onGUIEvent(self,controlID,valueName,value):
		print 'GUIEvent received: controldID='+controlID+' valueName='+valueName+' value='+value

		return 0 



	# key and mouse events; use this to add some user interaction to your scripts 
	def onKeyPressed(self,k):
		print 'onKeyPressed '+k

		return 0 

	def onKeyReleased(self,k):
		print 'onKeyReleased '+k

		return 0 

	def onMouseButtonLeft(self,x,y,pressed):
		print 'onMouseButtonLeft x='+str(x)+' y='+str(y)+' pressed='+str(pressed)

		return 0

	def onMouseButtonRight(self,x,y,pressed):
		print 'onMouseButtonRight x='+str(x)+' y='+str(y)+' pressed='+str(pressed)

		return 0

	def onMouseButtonMiddle(self,x,y,pressed):
		print 'onMouseButtonMiddle x='+str(x)+' y='+str(y)+' pressed='+str(pressed)

		return 0

	def onMouseWheel(self,x,y,delta):
		print 'onMouseButtonWheel x='+str(x)+' y='+str(y)+' delta='+str(delta)

		return 0




 
