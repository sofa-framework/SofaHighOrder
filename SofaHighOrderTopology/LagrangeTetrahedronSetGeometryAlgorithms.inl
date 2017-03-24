#ifndef SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_INL

#include "LagrangeTetrahedronSetGeometryAlgorithms.h"
#include "HighOrderTetrahedronSetGeometryAlgorithms.inl"
#include <sofa/core/visual/VisualParams.h>
#include "LagrangeTetrahedronCoefficients.h"
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/rmath.h>
#include <fstream>
namespace sofa
{

namespace component
{

namespace topology
{


const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};

 size_t moduleTVI(TetrahedronIndexVector tvi) {
	 return (tvi[0]+tvi[1]+tvi[2]+tvi[3]);
 }
 size_t factorialTVI(TetrahedronIndexVector tvi) {
	 return (lfactorial(tvi[0])*lfactorial(tvi[1])*lfactorial(tvi[2])*lfactorial(tvi[3]));
 }


template< class DataTypes>
 LagrangeTetrahedronSetGeometryAlgorithms< DataTypes >::LagrangeTetrahedronSetGeometryAlgorithms() : 
HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>()
    {
    }
 template< class DataTypes>
long int LagrangeTetrahedronSetGeometryAlgorithms< DataTypes >::sterlingVector(const TetrahedronIndexVector tbiIn,const TetrahedronIndexVector tbiIk) const
 {
	size_t i;
	long int val;
	val=1;
	for (i=0;i<4;++i) {
		val*=stirlingNumberArray[tbiIn[i]][tbiIk[i]];
	}
	return(val);
}
template< class DataTypes>
 void LagrangeTetrahedronSetGeometryAlgorithms< DataTypes >::init()
{
	HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::init();
	/// get the this->degree of the Lagrange tetrahedron
	this->degree=this->container->getDegree();
	/// store the tetrahedron vector  index for each tetrahedron
	this->tbiArray=this->container->getTetrahedronIndexArray();
	/// compute the Lagrange coefficient for each control point in a tetrahedron
	
	/// precompute the array of Stirling number of the (signed) first kind using a recursive function
	 for (size_t n=0;n<=this->degree;++n) {
		 helper::vector<Real> stArray(n+1);
		 stArray[0]=0.0;
		 for (size_t k=1;k<n;++k) {
			 stArray[k]=stirlingNumberArray[n-1][k-1]-(n-1)*stirlingNumberArray[n-1][k];
		 }
		 stArray[n]=(Real)1.0;
		 stirlingNumberArray.push_back(stArray);
//		 for (size_t k=0;k<=n;++k) 
//			 std::cerr<< "Stirling number ["<<n<<", "<<k<<"]="<<stirlingNumberArray[n][k]<<std::endl;
	 }
//	checkCoefficients();
}




template<class DataTypes>
typename DataTypes::Real LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunction(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
{
	size_t i,j;
	Real val=1;
	for (i=0;i<4;++i) {

		for (j=1;j<=tbi[i];++j) {
			val*=(this->degree*barycentricCoordinate[i]-j+1)/(Real)j;
		}

	}
	return(val);
}
 template<class DataTypes>
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Vec4 LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunctionDerivatives(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
 {
	 size_t i,j,k;
	 Vec4 dval(1,1,1,1);
	 for (i=0;i<4;++i) {

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
			 for(j=1;j<4;++j) {
				 dval[(i+j)%4]*=val2;
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
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunctionHessian(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
 {


	 size_t i,j,k;
	 Mat44 hessian;
	 /*
	 Vec4 value,dvalue;
	 /// compute the value and its derivative of each term of the shape function
	 for(i=0;i<4;++i) {
		 Real val=1;
		 Real dval=1;
		 Real valDenom=0;
  	     Real valZero=0;

		 for (j=1;j<tbi[i];++j) {
 			tmp=(this->degree*barycentricCoordinate[i]-j+1)/(Real)j;
			 val*=tmp; 
			  if (tmp==0) {
				 valZero=this->degree/(Real)(j);
			 } else {
				 dval*=tmp;
				 valDenom+=this->degree/(Real)(j*tmp);
			 }
		 }

		 value[i]=val;
		 if (valZero==0) 
			 dval=dval/valDenom;
		 else 
			 dval=dval*valZero;
		 dvalue[i]=dval;
	 }
	 // now compute the hessian for the off diagonal terms
	 for(i=0;i<4;++i) {
		 for(j=i+1;j<4;++j) {
			 hessian[i][j]=1;
			 for (k=0;k<4;++k) {
				 if ((k!=i)&&(k!=j)) {
					 hessian[i][j]*=value[k];
				 }
			 }
			 hessian[i][j]*=dvalue[i]*dvalue[j];
			  hessian[j][i]
		 }
	 }
			  Real  val=1;

			 hessian[i][j]=1.0;
		 }
	 }
	 
	 for (i=0;i<4;++i) {
		 Real  val=1;
		 Real tmp;
		 Real valDenom=0;
		 Real valZero=0;

		 for (j=1;j<tbi[i];++j) {
			 tmp=(this->degree*barycentricCoordinate[i]-j+1)/(Real)j;
			 
			 if (tmp==0) {
				 valZero=this->degree/(Real)(j);
			 } else {
				 val*=tmp;
				 valDenom+=this->degree/(Real)(j*tmp);
			 }
		 }
		 for(j=1;j<4;++j) {
			 dval[(i+j)%4]*=val;
		 }
		 if (valZero==0) 
			 dval[i]=val/valDenom;
		 else 
			 dval[i]=val*valZero;

	 } */
     return hessian;
 }
 
  template<class DataTypes>
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Real LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {

 	 TetrahedronIndexVector tbiu,tbiv;
	 size_t al,am,an,ai,aj,ak,jj,kk;
	 Real val2=6.0/((Real) factorialTVI(tbi1)* factorialTVI(tbi2));
	 Real val=0;
	 Real tmp;
	 size_t sumIndex=0;
	 for (al=0;al<=tbi1[0];++al) {
		 for (am=0;am<=tbi1[1];++am) {
			 for (an=0;an<=tbi1[2];++an) {
				 for (ai=0;ai<=tbi1[3];++ai) {
					 tbiu=TetrahedronIndexVector(al,am,an,ai);
					 for (aj=0;aj<=tbi2[0];++aj) {
						 for (ak=0;ak<=tbi2[1];++ak) {
							 for (jj=0;jj<=tbi2[2];++jj) {
								 for (kk=0;kk<=tbi2[3];++kk) {
									 tbiv=TetrahedronIndexVector(aj,ak,jj,kk); 
									 sumIndex= moduleTVI(tbiu)+moduleTVI(tbiv);
									 tmp=(Real)factorialTVI(tbiu+tbiv)/((Real)lfactorial(3+sumIndex));
									 tmp*=(Real)sterlingVector(tbi1,tbiu)*(Real)sterlingVector(tbi2,tbiv);  

									 val+=pow((Real)this->degree,(Real)sumIndex)*tmp;
								 }
							 }
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
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Real LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {
 size_t i=  (this->degree+1)*(this->degree+2)*(this->degree+3)/6-
		 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])*(this->degree+3-tbi1[0])/6+
		 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2-
		 (this->degree+1-tbi1[0]-tbi1[1])*(this->degree+2-tbi1[0]-tbi1[1])/2+ tbi1[2]; 
	 size_t j=(this->degree+1)*(this->degree+2)*(this->degree+3)/6-
		 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])*(this->degree+3-tbi2[0])/6+
		 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2-
		 (this->degree+1-tbi2[0]-tbi2[1])*(this->degree+2-tbi2[0]-tbi2[1])/2+ tbi2[2];
	 size_t ii=std::min(i,j);
     size_t jj=std::max(i,j);
	 size_t s=this->tbiArray.size();
	 switch(this->degree) 
	 {
	 case 2:
		 return(lagrangeTetrahedronDegree2AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 3:
		 return(lagrangeTetrahedronDegree3AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 4:
 		 return(lagrangeTetrahedronDegree4AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 5:
		 return(lagrangeTetrahedronDegree5AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 default:
		 {

			 Real val=computeAffineMassCoefficient(tbi1,tbi2);
			 return(val);
		 }
	 }
 }

 template<class DataTypes>
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Real LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,
	 const size_t indr[3])
 {
	 // compute the mass coefficient in closed form using an horrible expression
	 Real val=0;
	 unsigned int permutation[6][4]= {{0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {2, 1, 0, 3}};
	 unsigned int triplets[8][3] = {{1, 2, 3}, {1, 2, 4}, {1, 4, 4}, {4, 2, 3}, {4, 4, 3}, {4, 4, 4}, {1, 4, 3}, {4, 2, 4}};
	 int tripletSign[8] ={1, -1, 1, -1, 1, -1, -1, 1}; 
	 int sign3[6] = {1, -1, -1, 1, 1, -1};

	 size_t i,j,k;
	 TetrahedronIndexVector tbir[3],tbiu,tbiv,tbiw,tbix,tbiy,tbiz;
	 for (i=0;i<3;++i) {
		 tbir[i]=this->tbiArray[indr[i]];
	 }
	 unsigned int newIndex[3];

	 size_t aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,sumIndex;
	 Real denom=factorialTVI(tbi1)*factorialTVI(tbi2)*factorialTVI(tbir[0])*factorialTVI(tbir[1])*factorialTVI(tbir[2]);

	 val=0;
	 Real tmp;


	 for (i=0;i<6;++i) {

		 for(j=0;j<8;++j) {
			 for(k=0;k<3;++k) {
				 newIndex[k]=permutation[i][triplets[j][k]-1];
			 }
			 for (aa=0;aa<=tbi1[0];++aa) {
				 for (ab=0;ab<=tbi1[1];++ab) {
					 for (ac=0;ac<=tbi1[2];++ac) {
						 for (ad=0;ad<=tbi1[3];++ad) {
							 tbiu=TetrahedronIndexVector(aa,ab,ac,ad);
							 for (ae=0;ae<=tbi2[0];++ae) {
								 for (af=0;af<=tbi2[1];++af) {
									 for (ag=0;ag<=tbi2[2];++ag) {
										 for (ah=0;ah<=tbi2[3];++ah) {
											 tbiv=TetrahedronIndexVector(ae,af,ag,ah);
											 for (ai=0;ai<=tbir[0][0];++ai) {
												 for (aj=0;aj<=tbir[0][1];++aj) {
													 for (ak=0;ak<=tbir[0][2];++ak) {
														 for (al=0;al<=tbir[0][3];++al) {
															 tbiw=TetrahedronIndexVector(ai,aj,ak,al);
															 if (tbiw[newIndex[0]]>0) {
																 for (am=0;am<=tbir[1][0];++am) {
																	 for (an=0;an<=tbir[1][1];++an) {
																		 for (ao=0;ao<=tbir[1][2];++ao) {
																			 for (ap=0;ap<=tbir[1][3];++ap) {
																				 tbix=TetrahedronIndexVector(am,an,ao,ap);
																				 if (tbix[newIndex[1]]>0) {
																					 for (aq=0;aq<=tbir[2][0];++aq) {
																						 for (ar=0;ar<=tbir[2][1];++ar) {
																							 for (as=0;as<=tbir[2][2];++as) {
																								 for (at=0;at<=tbir[2][3];++at) {
																									 tbiy=TetrahedronIndexVector(aq,ar,as,at);
																									 if (tbiy[newIndex[2]]>0) {
																										 tmp=(Real)sterlingVector(tbi1,tbiu)*(Real)sterlingVector(tbi2,tbiv)*
																											 (Real)sterlingVector(tbir[0],tbiw)*(Real)sterlingVector(tbir[1],tbix)*
																											 (Real)sterlingVector(tbir[2],tbiy);

																										 sumIndex= moduleTVI(tbiu)+moduleTVI(tbiv)+moduleTVI(tbiw)+
																											 moduleTVI(tbix)+moduleTVI(tbiy);
																										 tmp*=tbiw[newIndex[0]]*tbix[newIndex[1]]*tbiy[newIndex[2]];
																										 tbiz=tbiu+tbiv+tbiw+tbix+tbiy;
																										 --tbiz[newIndex[0]];--tbiz[newIndex[1]];--tbiz[newIndex[2]];
																										 tmp*=factorialTVI(tbiz)/(Real)lfactorial(sumIndex);
																										 val+=tripletSign[j]*sign3[i]*pow((Real)this->degree,(Real)sumIndex)*tmp;
																										 //if (tmp!=0)  {
																										 // tmp*=tripletSign[j]*sign3[i]*pow((Real)this->degree,(Real)sumIndex);
																										 // std::cerr<< "tmp="<<tmp<<std::endl;
																										 //}
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
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Real LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::getRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,
	 const size_t indr[3])
 {
	 Real val=0;
	  if (this->degree<4) {
		 size_t k,r,s,t;

		 size_t Nd=(this->degree+1)*(this->degree+2)*(this->degree+3)/6;
		 size_t nbElem=Nd*(Nd-1)*(Nd-2)/6;
		 // read coefficient from array. but index is based on lexicographic order of TVI and must be sorted in ascendent order
		 size_t i=  (this->degree+1)*(this->degree+2)*(this->degree+3)/6-
			 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])*(this->degree+3-tbi1[0])/6+
			 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2-
			 (this->degree+1-tbi1[0]-tbi1[1])*(this->degree+2-tbi1[0]-tbi1[1])/2+ tbi1[2]; 
		 size_t j=(this->degree+1)*(this->degree+2)*(this->degree+3)/6-
			 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])*(this->degree+3-tbi2[0])/6+
			 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2-
			 (this->degree+1-tbi2[0]-tbi2[1])*(this->degree+2-tbi2[0]-tbi2[1])/2+ tbi2[2];


		 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
		 std::vector<size_t> sub(3);
		 for (k=0;k<3;++k) 
			 sub[k]=this->container->getLexicographicIndex(indr[k]);
		 std::sort(sub.begin(),sub.end());
		 r=sub[0];s=sub[1];t=sub[2];

		 switch(this->degree) {
		 case 2:
			 val=lagrangeTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
			 break;
		 case 3:
			 val=lagrangeTetrahedronDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
			 break;
		 }
	  } else {
		 val=computeRegularMassCoefficient(tbi1,tbi2,indr);
	  }
	return val;
 }
 
 template<class DataTypes>
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::computeAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {
	 Mat44 mat;
	 TetrahedronIndexVector tbiu,tbiv,tbiz;
	 size_t al,am,an,ai,aj,ak,jj,kk,m,n;
	 Real val2=6.0/((Real) factorialTVI(tbi1)* factorialTVI(tbi2));
	 Real val=0;
	 Real tmp;
	 size_t sumIndex=0;
	 for (al=0;al<=tbi1[0];++al) {
		 for (am=0;am<=tbi1[1];++am) {
			 for (an=0;an<=tbi1[2];++an) {
				 for (ai=0;ai<=tbi1[3];++ai) {
					 tbiu=TetrahedronIndexVector(al,am,an,ai);
					 for (aj=0;aj<=tbi2[0];++aj) {
						 for (ak=0;ak<=tbi2[1];++ak) {
							 for (jj=0;jj<=tbi2[2];++jj) {
								 for (kk=0;kk<=tbi2[3];++kk) {
									 tbiv=TetrahedronIndexVector(aj,ak,jj,kk);

									 for(m=0;m<4;++m) {
										 for(n=0;n<4;++n) {
											 if ((tbiu[m]>0) && (tbiv[n]>0)) {
												 sumIndex= moduleTVI(tbiu)+moduleTVI(tbiv);
												 tbiz=tbiu+tbiv;
												 --tbiz[m];--tbiz[n];
												 tmp=(Real)factorialTVI(tbiz)/((Real)lfactorial(1+sumIndex));
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
		 }
	 }
	 mat*=val2;


	 return(mat);
 }
  template<class DataTypes>
 typename LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {
	 Mat44 val;
	 size_t i=  (this->degree+1)*(this->degree+2)*(this->degree+3)/6-
		 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])*(this->degree+3-tbi1[0])/6+
		 (this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2-
		 (this->degree+1-tbi1[0]-tbi1[1])*(this->degree+2-tbi1[0]-tbi1[1])/2+ tbi1[2]; 
	 size_t j=(this->degree+1)*(this->degree+2)*(this->degree+3)/6-
		 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])*(this->degree+3-tbi2[0])/6+
		 (this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2-
		 (this->degree+1-tbi2[0]-tbi2[1])*(this->degree+2-tbi2[0]-tbi2[1])/2+ tbi2[2];
	 size_t ii=std::min(i,j);
	 size_t jj=std::max(i,j);
	 size_t k;
	 size_t s=this->tbiArray.size();
	 size_t index=(ii*(2*s-ii+1)/2+(jj-ii))*16;
	 bool transpose=false;
	 if (ii==j)
		 transpose=true;
	 switch(this->degree) 
	 {
	 case 2:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=lagrangeTetrahedronDegree2AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 3:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=lagrangeTetrahedronDegree3AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 4:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=lagrangeTetrahedronDegree4AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 5:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=lagrangeTetrahedronDegree5AffineStiffnessCoefficient[index+k];
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
 real oneNorm(const sofa::defaulttype::Mat<4,4, real>& A)
{
     real norm = 0.0;
    for (int i=0; i<4; i++)
    {
        real columnAbsSum = helper::rabs(A(0,i)) + helper::rabs(A(1,i)) + helper::rabs(A(2,i))+ helper::rabs(A(3,i));
        if (columnAbsSum > norm)
            norm = columnAbsSum;
    }
    return norm;
}
  template<class DataTypes>
void  LagrangeTetrahedronSetGeometryAlgorithms<DataTypes>::checkCoefficients()
 { 
	 // check value of shape function
	 {
		 size_t iii=4;
		 size_t i;
		 Real val;
		 Vec4 bc=((Vec4)this->tbiArray[iii]/(Real)this->degree);
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
		 size_t iii=4;
		 size_t i;
		 Vec4 val;
		 Vec4 bc=((Vec4)this->tbiArray[iii]/(Real)this->degree);
		 for (i=0;i<this->tbiArray.size();++i) {
			 val=computeShapeFunctionDerivatives(this->tbiArray[this->container->getHierarchicalIndex(i)],bc);
//			 if (i==iii)
//				 assert((val[0]==(0.0))&& (val[1]==(0.0))&& (val[2]==(0.0))&& (val[3]==(0.0)));
		 }
	 }

	 // number of control points
	 size_t Nd=(this->degree+1)*(this->degree+2)*(this->degree+3)/6;
	 // number of unique entries in the symmetric mass matrix
	 size_t nbMassEntries=Nd*(Nd+1)/2;
	 if (this->degree<6) {
		 // check mass coefficients
		 size_t i,j,k,l,m;
		 size_t s=this->tbiArray.size();
		 size_t nbStiffnessEntries=nbMassEntries*16;
		 for (i=0;i<this->tbiArray.size();++i) {

			 size_t ai=this->container->getHierarchicalIndex(i);
			 for (j=i;j<this->tbiArray.size();++j) {
				 size_t index=(i*(2*s-i+1)/2+(j-i))*16;
				 size_t aj=this->container->getHierarchicalIndex(j);
				 Mat44 val=computeAffineStiffnessCoefficientMatrix(this->tbiArray[ai],this->tbiArray[aj]);
				 Mat44 val2;
				 switch(this->degree) {
				 case 2:
					 assert(nbStiffnessEntries==arraysize(lagrangeTetrahedronDegree2AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=lagrangeTetrahedronDegree2AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 3:
					 assert(nbStiffnessEntries==arraysize(lagrangeTetrahedronDegree3AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=lagrangeTetrahedronDegree3AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 4:
					 assert(nbStiffnessEntries==arraysize(lagrangeTetrahedronDegree4AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=lagrangeTetrahedronDegree4AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 5:
					 assert(nbStiffnessEntries==arraysize(lagrangeTetrahedronDegree5AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=lagrangeTetrahedronDegree5AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
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
				 TetrahedronIndexVector tbi1,tbi2;
				 tbi1=this->tbiArray[aai];
				 tbi2=this->tbiArray[aaj];
				 Real val=computeAffineMassCoefficient(tbi1,tbi2);

				 switch(this->degree) {
				 case 2:
					 assert(nbMassEntries==arraysize(lagrangeTetrahedronDegree2AffineMassCoefficient)); 
					 assert(fabs(lagrangeTetrahedronDegree2AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 3:
					 assert(nbMassEntries==arraysize(lagrangeTetrahedronDegree3AffineMassCoefficient)); 
					 assert(fabs(lagrangeTetrahedronDegree3AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 4:
					 assert(nbMassEntries==arraysize(lagrangeTetrahedronDegree4AffineMassCoefficient)); 
					 assert(fabs(lagrangeTetrahedronDegree4AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 5:
					 assert(nbMassEntries==arraysize(lagrangeTetrahedronDegree5AffineMassCoefficient)); 
					 assert(fabs(lagrangeTetrahedronDegree5AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 }

			 }
		 } 	 
	 }
	 if (this->degree<3) {
		 
		 size_t i,j,r,s,t;
		 size_t Ndd=(this->degree+1)*(this->degree+2)*(this->degree)/6;
		 size_t nbElem=Nd*(Nd-1)*(Nd-2)/6;
		 size_t sss=this->tbiArray.size();
		 size_t aSize=arraysize(lagrangeTetrahedronDegree2RegularMassCoefficient);
		 assert(aSize==(nbElem*nbMassEntries));

		 for (i=0;i<this->tbiArray.size();++i) {

			 size_t ai=this->container->getHierarchicalIndex(i);
			 for (j=i;j<this->tbiArray.size();++j) {
				  size_t aj=this->container->getHierarchicalIndex(j);
					
				
				 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
				 size_t sub[3];
				 for (r=0;r<this->tbiArray.size();++r) {
					 size_t ar=this->container->getHierarchicalIndex(r);
					 for (s=r+1;s<this->tbiArray.size();++s) {
						 size_t as=this->container->getHierarchicalIndex(s);
						 for (t=s+1;t<this->tbiArray.size();++t) {
							size_t at=this->container->getHierarchicalIndex(t);
							 sub[0]=ar;sub[1]=as;sub[2]=at;
							 Real val=computeRegularMassCoefficient(this->tbiArray[ai],this->tbiArray[aj],sub);
		
							 Real val2=lagrangeTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
							  assert(fabs(lagrangeTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)]-val)<1e-8);
						 }
					 }
				 }
			 }
		 } 
	 }
	 //else {
		// size_t i,j,r,s,t;
		// size_t Ndd=(this->degree+1)*(this->degree+2)*(this->degree)/6;
		// size_t nbElem=Nd*(Nd-1)*(Nd-2)/6;
		// size_t sss=this->tbiArray.size();

 	//	 std::vector<Real> coeffArray(nbElem*nbMassEntries);

		// for (i=0;i<this->tbiArray.size();++i) {

		//	 size_t ai=this->container->getHierarchicalIndex(i);
		//	 for (j=i;j<this->tbiArray.size();++j) {
		//		  size_t aj=this->container->getHierarchicalIndex(j);
		//			
		//		
		//		 size_t massIndex=(i*(2*Nd-i+1)/2+(j-i))*nbElem;
		//		 size_t sub[3];
		//		 for (r=0;r<this->tbiArray.size();++r) {
		//			 size_t ar=this->container->getHierarchicalIndex(r);
		//			 for (s=r+1;s<this->tbiArray.size();++s) {
		//				 size_t as=this->container->getHierarchicalIndex(s);
		//				 for (t=s+1;t<this->tbiArray.size();++t) {
		//					 size_t at=this->container->getHierarchicalIndex(t);
		//					 sub[0]=ar;sub[1]=as;sub[2]=at;
		//					 Real val=computeRegularMassCoefficient(this->tbiArray[ai],this->tbiArray[aj],sub);
		//					 coeffArray[massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)]=val;

		//				 }
		//			 }
		//		 }
		//	 }
		// } 
		// // write array in file
		// std::ostringstream st;
		// st<< "LagrangeTetrahedronDegree"<<(int)this->degree<<"RegMass.h";
		// std::string fname=st.str();

		// std::ofstream myfile("toto.h");
		// myfile.precision(16);
		// myfile <<"double lagrangeTriangleDegree"<<(int)this->degree<<"RegularMassCoefficient[]={";
		// myfile<<coeffArray[0];
		// for (i=1;i<coeffArray.size();++i) {
		//	 myfile<<","<<coeffArray[i];
		// }
		// myfile<<"};"<<std::endl;
		// myfile.close();
	 //}
}



} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_TETEAHEDRONSETGEOMETRYALGORITHMS_INL
