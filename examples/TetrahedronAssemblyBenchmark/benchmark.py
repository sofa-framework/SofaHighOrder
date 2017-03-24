import os,sys

FILE="benchmark.params"
CMDLINE="runSofa benchmarkTetrahedronElasticity.xml -g batch -n 0"

degree=2
resolution=2
nbTrees=6
componentType="elasticityTetra"
integrationType="analytical"
numericalIntegrationOrder=4


resolutionArrayTetra=[24,12,8,6,5]
resolutionArrayTrian=[72,36,25,18,15]
integrationOrderArray=[1,4,6,8,10]
integrationOrderArrayLumped=[1,2,3,4,5]



componentTypeList=['massTetra','massTrian','massTetraLumped','massTrianLumped']
integrationTypeList=['numerical','analytical','exact']		
for componentType in componentTypeList:
	print " componentType=",componentType
	for integrationType in integrationTypeList:
		print " integrationType=",integrationType
		if (integrationTypeList.index(integrationType)<2):
			rr=range(1,6)
		else:
			rr=range(2,4)
			if (componentTypeList.index(componentType)%2==1):
				rr=range(2,6)
		for degree in rr:
			print " degree=",degree
			sys.stdout.flush()
			print "index=",componentTypeList.index(componentType)
			if (componentTypeList.index(componentType)>1):
				order=integrationOrderArrayLumped[degree-1]
			else:
				order=integrationOrderArray[degree-1]
			if (componentTypeList.index(componentType)%2==0):
				resolution=resolutionArrayTetra[degree-1]
			else:
				resolution=resolutionArrayTrian[degree-1]
			print " resolution=",resolution				
			f=open(FILE,"w")
			line=str(degree)+" "+str(nbTrees)+" "+str(resolution)+" "+componentType+" "+integrationType+" "+str(order)
			f.write(line)
			f.close()
			os.system(CMDLINE)	
			
integrationOrderArrayTrian=[1,2,4,6,8]
integrationOrderArrayTetra=[1,2,4,6,8]			
	
componentTypeList=['elasticityTetra','diffusionTetra','diffusionTetraAniso','elasticityTrian','diffusionTrian','diffusionTrianAniso']
integrationTypeList=['numerical','standard','analytical']
for componentType in componentTypeList:
	print " componentType=",componentType
	for integrationType in integrationTypeList:
		print " integrationType=",integrationType
		for degree in range(1,6):
			print " degree=",degree
			sys.stdout.flush()
			print "index=",componentTypeList.index(componentType)
			if (componentTypeList.index(componentType)>2):
				resolution=resolutionArrayTrian[degree-1]
				order=integrationOrderArrayTrian[degree-1]
			else:
				resolution=resolutionArrayTetra[degree-1]
				order=integrationOrderArrayTetra[degree-1]
			print " resolution=",resolution				
			f=open(FILE,"w")
			line=str(degree)+" "+str(nbTrees)+" "+str(resolution)+" "+componentType+" "+integrationType+" "+str(order)
			f.write(line)
			f.close()
			os.system(CMDLINE)
			
