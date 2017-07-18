#ifndef SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_INL

#include "LagrangeTriangleSetGeometryAlgorithms.h"
#include "HighOrderTriangleSetGeometryAlgorithms.inl"
#include <sofa/core/visual/VisualParams.h>
#include "LagrangeTriangleCoefficients.h"
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/rmath.h>
#include <fstream>
namespace sofa
{

namespace component
{

namespace topology
{






 size_t moduleTVI3(TriangleIndexVector tvi) {
	 return (tvi[0]+tvi[1]+tvi[2]);
}

template< class DataTypes>
 LagrangeTriangleSetGeometryAlgorithms< DataTypes >::LagrangeTriangleSetGeometryAlgorithms() : 
HighOrderTriangleSetGeometryAlgorithms<DataTypes>()
    {
    }
 template< class DataTypes>
long int LagrangeTriangleSetGeometryAlgorithms< DataTypes >::sterlingVector(const TriangleIndexVector tbiIn,const TriangleIndexVector tbiIk) const
 {
	size_t i;
	long int val;
	val=1;
	for (i=0;i<3;++i) {
		val*=stirlingNumberArray[tbiIn[i]][tbiIk[i]];
	}
	return(val);
}

template< class DataTypes>
 void LagrangeTriangleSetGeometryAlgorithms< DataTypes >::init()
{
	HighOrderTriangleSetGeometryAlgorithms<DataTypes>::init();
	/// get the this->degree of the Lagrange triangle
	this->degree=this->container->getDegree();
	/// store the triangle vector  index for each triangle
	this->tbiArray=this->container->getTriangleIndexArray();
	/// compute the Lagrange coefficient for each control point in a triangle
	
	/// precompute the array of Stirling number of the (signed) first kind using a recursive function
	 for (size_t n=0;n<=this->degree;++n) {
		 helper::vector<Real> stArray(n+1);
		 stArray[0]=0.0;
		 for (size_t k=1;k<n;++k) {
			 stArray[k]=stirlingNumberArray[n-1][k-1]-(n-1)*stirlingNumberArray[n-1][k];
		 }
		 stArray[n]=(Real)1.0;
		 stirlingNumberArray.push_back(stArray);

	 }
//	checkCoefficients();
}




template<class DataTypes>
typename DataTypes::Real LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunction(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
{
	size_t i,j;
	Real val=1;
	for (i=0;i<3;++i) {

		for (j=1;j<=tbi[i];++j) {
			val*=(this->degree*barycentricCoordinate[i]-j+1)/(Real)j;
		}

	}
	return(val);
}

 template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Vec3 LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionDerivatives(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
 {
	 size_t i,j,k;
	 Vec3 dval(1,1,1);
	 for (i=0;i<3;++i) {

		 if (tbi[i]==0) {
			  dval[i]=0;
		 } else {
			 Real val=1;
			 Real val2=1;
			 Real tmp;
			 Real valDenom=0;
			 Real valZero=0;
			 for (j=1;j<=tbi[i];++j) {
				 tmp=(this->degree*barycentricCoordinate[i]-j+1)/(Real)j;
				 val2*=tmp;
				 if (tmp==0) {
					 valZero=(Real)this->degree/(Real)(j);

				 } else {
					 val*=tmp;
					 valDenom+=(Real)this->degree/(Real)(j*tmp);
				 }
			 }
			 for(j=1;j<3;++j) {
				 dval[(i+j)%3]*=val2;
			 }
			 if (valZero==0) {
				 dval[i]*=val*valDenom;
			 }
			 else 
				 dval[i]*=val*valZero;
		 }
	 }
     return dval;
 }



 template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Mat33 LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionHessian(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
 {


	 size_t i,j,k;
	 Mat33 hessian;

     return hessian;
 }

  template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Real LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {

 	 TriangleIndexVector tbiu,tbiv;
	 size_t al,am,an,aj,ak,jj;
	 Real val2=2.0/((Real) factorialTVI(tbi1)* factorialTVI(tbi2));
	 Real val=0;
	 Real tmp;
	 size_t sumIndex=0;
	 for (al=0;al<=tbi1[0];++al) {
		 for (am=0;am<=tbi1[1];++am) {
			 for (an=0;an<=tbi1[2];++an) {

				 tbiu=TriangleIndexVector(al,am,an);
				 for (aj=0;aj<=tbi2[0];++aj) {
					 for (ak=0;ak<=tbi2[1];++ak) {
						 for (jj=0;jj<=tbi2[2];++jj) {

							 tbiv=TriangleIndexVector(aj,ak,jj); 
							 sumIndex= moduleTVI3(tbiu)+moduleTVI3(tbiv);
							 tmp=(Real)factorialTVI(tbiu+tbiv)/((Real)lfactorial(2+sumIndex));
							 tmp*=(Real)sterlingVector(tbi1,tbiu)*(Real)sterlingVector(tbi2,tbiv);  

							 val+=pow((Real)this->degree,(Real)sumIndex)*tmp;
						 }

					 }
				 }
			 }
		 }
	 }
	 val*=val2;
	 return(val);
 }
 template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Real LagrangeTriangleSetGeometryAlgorithms<DataTypes>::getAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	 size_t i=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2+  tbi1[1]; 
	 size_t j=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2+  tbi2[1]; 
	 size_t ii=std::min(i,j);
     size_t jj=std::max(i,j);
	 size_t s=this->tbiArray.size();
	 switch(this->degree) 
	 {
	 case 2:
		 return(lagrangeTriangleDegree2AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 3:
		 return(lagrangeTriangleDegree3AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 4:
 		 return(lagrangeTriangleDegree4AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 5:
		 return(lagrangeTriangleDegree5AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 default:
		 {

			 Real val=computeAffineMassCoefficient(tbi1,tbi2);
			 return(val);
		 }
	 }
 }

 template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Real LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,
	 const size_t indr[2])
 {
	 // compute the mass coefficient in closed form using an horrible expression
	 Real val=0;
	 unsigned int permutation[2][3]= {{0, 1, 2}, {1, 0, 2}};
	 unsigned int triplets[4][2] = {{1, 2}, {1, 3}, {3, 2}, {3, 3}};
	 int tripletSign[4] ={1, -1, -1, 1}; 
	 int sign3[2] = {1, -1};

	 size_t i,j,k;
	 TriangleIndexVector tbir[3],tbiu,tbiv,tbiw,tbix,tbiy,tbiz;
	 for (i=0;i<2;++i) {
		 tbir[i]=this->tbiArray[indr[i]];
	 }
	 unsigned int newIndex[3];

	 size_t aa,ab,ac,ae,af,ag,ai,aj,ak,am,an,ao,sumIndex;
	 Real denom=factorialTVI(tbi1)*factorialTVI(tbi2)*factorialTVI(tbir[0])*factorialTVI(tbir[1]);

	 val=0;
	 Real tmp;


	 for (i=0;i<2;++i) {

		 for(j=0;j<4;++j) {
			 for(k=0;k<2;++k) {
				 newIndex[k]=permutation[i][triplets[j][k]-1];
			 }
			 for (aa=0;aa<=tbi1[0];++aa) {
				 for (ab=0;ab<=tbi1[1];++ab) {
					 for (ac=0;ac<=tbi1[2];++ac) {

						 tbiu=TriangleIndexVector(aa,ab,ac);
						 for (ae=0;ae<=tbi2[0];++ae) {
							 for (af=0;af<=tbi2[1];++af) {
								 for (ag=0;ag<=tbi2[2];++ag) {

									 tbiv=TriangleIndexVector(ae,af,ag);
									 for (ai=0;ai<=tbir[0][0];++ai) {
										 for (aj=0;aj<=tbir[0][1];++aj) {
											 for (ak=0;ak<=tbir[0][2];++ak) {

												 tbiw=TriangleIndexVector(ai,aj,ak);
												 if (tbiw[newIndex[0]]>0) {
													 for (am=0;am<=tbir[1][0];++am) {
														 for (an=0;an<=tbir[1][1];++an) {
															 for (ao=0;ao<=tbir[1][2];++ao) {

																 tbix=TriangleIndexVector(am,an,ao);
																 if (tbix[newIndex[1]]>0) {

																	 tmp=(Real)sterlingVector(tbi1,tbiu)*(Real)sterlingVector(tbi2,tbiv)*
																		 (Real)sterlingVector(tbir[0],tbiw)*(Real)sterlingVector(tbir[1],tbix);


																	 sumIndex= moduleTVI3(tbiu)+moduleTVI3(tbiv)+moduleTVI3(tbiw)+
																		 moduleTVI3(tbix);
																	 tmp*=tbiw[newIndex[0]]*tbix[newIndex[1]];
																	 tbiz=tbiu+tbiv+tbiw+tbix;
																	 --tbiz[newIndex[0]];--tbiz[newIndex[1]];
																	 tmp*=factorialTVI(tbiz)/(Real)lfactorial(sumIndex);
																	 val+=tripletSign[j]*sign3[i]*pow((Real)this->degree,(Real)sumIndex)*tmp;

																 }
															 }
														 }
													 }
												 }
											 }
										 }

									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }

	 }
 
	 val/=denom;
	 return(val);
 }

template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Real LagrangeTriangleSetGeometryAlgorithms<DataTypes>::getRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,
	 const size_t indr[2])
 {
	 Real val=0;
	  if (this->degree<6) {
		 size_t k,r,s;

		 size_t Nd=(this->degree+1)*(this->degree+2)/2;
		 size_t nbElem=Nd*(Nd-1)/2;
		 // read coefficient from array. but index is based on lexicographic order of TVI and must be sorted in ascendent order
		 size_t i=(this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2+  tbi1[1]; 
		 size_t j=(this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2+  tbi2[1]; 


		 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
		 std::vector<size_t> sub(2);
		 for (k=0;k<2;++k) 
			 sub[k]=this->container->getLexicographicIndex(indr[k]);
		 std::sort(sub.begin(),sub.end());
		 r=sub[0];s=sub[1];
		 switch(this->degree) {
		 case 2:
			 val=lagrangeTriangleDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];

			 break;
		 case 3:
			 val=lagrangeTriangleDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 case 4:
			 val=lagrangeTriangleDegree4RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 case 5:
			 val=lagrangeTriangleDegree5RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 }

	  } else {
		 val=computeRegularMassCoefficient(tbi1,tbi2,indr);
	  }
	return val;
 }
 
 template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Mat33 LagrangeTriangleSetGeometryAlgorithms<DataTypes>::computeAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	 Mat33 mat;
	 TriangleIndexVector tbiu,tbiv,tbiz;
	 size_t al,am,an,aj,ak,jj,m,n;
	 Real val2=2.0/((Real) factorialTVI(tbi1)* factorialTVI(tbi2));
	 Real val=0;
	 Real tmp;
	 size_t sumIndex=0;
	 for (al=0;al<=tbi1[0];++al) {
		 for (am=0;am<=tbi1[1];++am) {
			 for (an=0;an<=tbi1[2];++an) {

				 tbiu=TriangleIndexVector(al,am,an);
				 for (aj=0;aj<=tbi2[0];++aj) {
					 for (ak=0;ak<=tbi2[1];++ak) {
						 for (jj=0;jj<=tbi2[2];++jj) {

							 tbiv=TriangleIndexVector(aj,ak,jj);

							 for(m=0;m<3;++m) {
								 for(n=0;n<3;++n) {
									 if ((tbiu[m]>0) && (tbiv[n]>0)) {
										 sumIndex= moduleTVI3(tbiu)+moduleTVI3(tbiv);
										 tbiz=tbiu+tbiv;
										 --tbiz[m];--tbiz[n];
										 tmp=(Real)factorialTVI(tbiz)/((Real)lfactorial(sumIndex));
										 tmp*=tbiu[m]*tbiv[n];
										 tmp*=(Real)sterlingVector(tbi1,tbiu)*(Real)sterlingVector(tbi2,tbiv);  

										 mat[m][n]+=pow((Real)this->degree,(Real)sumIndex)*tmp;
									 }
								 }
							 }
						 }
					 }
				 }
			 }

		 }
	 }
	 mat*=val2;


	 return(mat);
 }
  template<class DataTypes>
 typename LagrangeTriangleSetGeometryAlgorithms<DataTypes>::Mat33 LagrangeTriangleSetGeometryAlgorithms<DataTypes>::getAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	 Mat33 val;
	 size_t i=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2+  tbi1[1]; 
	 size_t j=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2+  tbi2[1]; 
	 size_t ii=std::min(i,j);
	 size_t jj=std::max(i,j);
	 size_t k;
	 size_t s=this->tbiArray.size();
	 size_t index=(ii*(2*s-ii+1)/2+(jj-ii))*9;
	 bool transpose=false;
	 if (ii==j)
		 transpose=true;
	 switch(this->degree) 
	 {
	 case 2:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=lagrangeTriangleDegree2AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 3:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=lagrangeTriangleDegree3AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 4:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=lagrangeTriangleDegree4AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 5:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=lagrangeTriangleDegree5AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 default:
		 val=computeAffineStiffnessCoefficientMatrix(tbi1,tbi2);
	 }	
	 if (transpose)
		 val.transpose();
	 return(val);
 }

  // good
template <typename T,unsigned S>
inline unsigned arraysize(const T (&v)[S]) { return S; }

template <typename  real>
 real oneNorm3(const sofa::defaulttype::Mat<3,3, real>& A)
{
     real norm = 0.0;
    for (int i=0; i<3; i++)
    {
        real columnAbsSum = helper::rabs(A(0,i)) + helper::rabs(A(1,i)) + helper::rabs(A(2,i));
        if (columnAbsSum > norm)
            norm = columnAbsSum;
    }
    return norm;
}
  template<class DataTypes>
void  LagrangeTriangleSetGeometryAlgorithms<DataTypes>::checkCoefficients()
 { 
	 // check value of shape function
	 {
		 size_t iii=3;
		 size_t i;
		 Real val;
		 Vec3 bc=((Vec3)this->tbiArray[iii]/(Real)this->degree);
		 for (i=0;i<this->tbiArray.size();++i) {
			 val=computeShapeFunction(this->tbiArray[i],bc);
			 if (i==iii)
				 assert(val==(1.0));
			 else 
				 assert(val==0.0);
		 }
	 }
	 // check value of shape function derivative
	 {
		 size_t iii=3;
		 size_t i;
		 Vec3 val;
		 Vec3 bc=((Vec3)this->tbiArray[iii]/(Real)this->degree);
		 for (i=0;i<this->tbiArray.size();++i) {
			 val=computeShapeFunctionDerivatives(this->tbiArray[this->container->getHierarchicalIndex(i)],bc);
//			 if (i==iii)
//				 assert((val[0]==(0.0))&& (val[1]==(0.0))&& (val[2]==(0.0))&& (val[3]==(0.0)));
		 }
	 }

	 // number of control points
	 size_t Nd=(this->degree+1)*(this->degree+2)/2;
	 // number of unique entries in the symmetric mass matrix
	 size_t nbMassEntries=Nd*(Nd+1)/2;
	 if (this->degree<6) {
		 // check mass coefficients
		 size_t i,j,k,l,m;
		 size_t s=this->tbiArray.size();
		 size_t nbStiffnessEntries=nbMassEntries*9;
		
		 for (i=0;i<this->tbiArray.size();++i) {

			 size_t ai=this->container->getHierarchicalIndex(i);
			 for (j=i;j<this->tbiArray.size();++j) {
				 size_t index=(i*(2*s-i+1)/2+(j-i))*9;
				 size_t aj=this->container->getHierarchicalIndex(j);
				 Mat33 val=computeAffineStiffnessCoefficientMatrix(this->tbiArray[ai],this->tbiArray[aj]);
				 Mat33 val2;
				 switch(this->degree) {
				 case 2:
					 assert(nbStiffnessEntries==arraysize(lagrangeTriangleDegree2AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=lagrangeTriangleDegree2AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm3(val-val2))<1e-8);
					 break;
				 case 3:
					 assert(nbStiffnessEntries==arraysize(lagrangeTriangleDegree3AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=lagrangeTriangleDegree3AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm3(val-val2))<1e-8);
					 break;
				 case 4:
					 assert(nbStiffnessEntries==arraysize(lagrangeTriangleDegree4AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=lagrangeTriangleDegree4AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm3(val-val2))<1e-8);
					 break;
				 case 5:
					 assert(nbStiffnessEntries==arraysize(lagrangeTriangleDegree5AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=lagrangeTriangleDegree5AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm3(val-val2))<1e-8);
					 break;
				 }

			 }
		 } 	 


	 }
	 if (this->degree<6) {
		 // check mass coefficients
		 size_t i,j;
		 size_t s=this->tbiArray.size();

		 for (i=0;i<this->tbiArray.size();++i) {

			 size_t aai=this->container->getHierarchicalIndex(i);
			 for (j=i;j<this->tbiArray.size();++j) {
				 size_t aaj=this->container->getHierarchicalIndex(j);
				 TriangleIndexVector tbi1,tbi2;
				 tbi1=this->tbiArray[aai];
				 tbi2=this->tbiArray[aaj];
				 Real val=computeAffineMassCoefficient(tbi1,tbi2);

				 switch(this->degree) {
				 case 2:
					 assert(nbMassEntries==arraysize(lagrangeTriangleDegree2AffineMassCoefficient)); 
					 assert(fabs(lagrangeTriangleDegree2AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 3:
					 assert(nbMassEntries==arraysize(lagrangeTriangleDegree3AffineMassCoefficient)); 
					 assert(fabs(lagrangeTriangleDegree3AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 4:
					 assert(nbMassEntries==arraysize(lagrangeTriangleDegree4AffineMassCoefficient)); 
					 assert(fabs(lagrangeTriangleDegree4AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 5:
					 assert(nbMassEntries==arraysize(lagrangeTriangleDegree5AffineMassCoefficient)); 
					 assert(fabs(lagrangeTriangleDegree5AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 }

			 }
		 } 	 
	 }
	 if (this->degree<6) { 
		 
		 size_t i,j,r,s;
		 size_t Ndd=(this->degree+1)*(this->degree+2)/2;
		 size_t nbElem=Nd*(Nd-1)/2;
		 size_t sss=this->tbiArray.size();
		 size_t aSize;
		 switch(this->degree) {
		 case 2:
			  aSize=arraysize(lagrangeTriangleDegree2RegularMassCoefficient);
			  break;
		 case 3:
			  aSize=arraysize(lagrangeTriangleDegree3RegularMassCoefficient);
			  break;
		 case 4:
			  aSize=arraysize(lagrangeTriangleDegree4RegularMassCoefficient);
			  break;
		 case 5:
			  aSize=arraysize(lagrangeTriangleDegree5RegularMassCoefficient);
			  break;
		 }
		 assert(aSize==(nbElem*nbMassEntries));

		 for (i=0;i<this->tbiArray.size();++i) {

			 size_t ai=this->container->getHierarchicalIndex(i);
			 for (j=i;j<this->tbiArray.size();++j) {
				 size_t aj=this->container->getHierarchicalIndex(j);


				 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
				 size_t sub[2];
				 for (r=0;r<this->tbiArray.size();++r) {
					 size_t ar=this->container->getHierarchicalIndex(r);
					 for (s=r+1;s<this->tbiArray.size();++s) {
						 size_t as=this->container->getHierarchicalIndex(s);

						 sub[0]=ar;sub[1]=as;
						 Real val=computeRegularMassCoefficient(this->tbiArray[ai],this->tbiArray[aj],sub);
						 Real val2;
						 switch(this->degree) {
						 case 2:
							 val2=lagrangeTriangleDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];

							 break;
						 case 3:
							 val2=lagrangeTriangleDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 case 4:
							 val2=lagrangeTriangleDegree4RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 case 5:
							 val2=lagrangeTriangleDegree5RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 }
						 assert(fabs(val2-val)<1e-8);
					 }
				 }
			 }
		 } 
	 } 
	// else {
	//	 std::cerr << "now creating the missing coefficients"<<std::endl;
	//	  size_t i,j,r,s;
	//	 size_t Ndd=(this->degree+1)*(this->degree+2)/2;
	//	 size_t nbElem=Nd*(Nd-1)/2;
	//	 size_t sss=this->tbiArray.size();
	//	 size_t aSize;

	//	 std::vector<Real> coeffArray(nbElem*nbMassEntries);


	//	 for (i=0;i<this->tbiArray.size();++i) {

	//		 size_t ai=this->container->getHierarchicalIndex(i);
	//		 for (j=i;j<this->tbiArray.size();++j) {
	//			 size_t aj=this->container->getHierarchicalIndex(j);


	//			 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
	//			 size_t sub[2];
	//			 for (r=0;r<this->tbiArray.size();++r) {
	//				 size_t ar=this->container->getHierarchicalIndex(r);
	//				 for (s=r+1;s<this->tbiArray.size();++s) {
	//					 size_t as=this->container->getHierarchicalIndex(s);

	//					 sub[0]=ar;sub[1]=as;
	//					 Real val=computeRegularMassCoefficient(this->tbiArray[ai],this->tbiArray[aj],sub);
	//					 coeffArray[massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)]=val;

	//				 }
	//			 }
	//		 }
	//	 } 
	//// write array in file
	//	 std::ostringstream st;
	//	 st<< "LagrangeTriangleDegree"<<this->degree<<"RegMass.h";
	//	 std::string fname=st.str();

	//	 std::ofstream myfile("toto,h");
	//	 myfile.precision(16);
	//	 myfile <<"double lagrangeTriangleDegree"<<this->degree<<"RegularMassCoefficient[]={";
	//	 myfile<<coeffArray[0];
	//	 for (i=1;i<coeffArray.size();++i) {
	//		 myfile<<","<<coeffArray[i];
	//	 }
	//	 myfile<<"};"<<std::endl;
	//	 myfile.close();

	// }
}


} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_TETEAHEDRONSETGEOMETRYALGORITHMS_INL
