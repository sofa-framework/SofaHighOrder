#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_INL

#include "BezierTetrahedronSetGeometryAlgorithms.h"
#include "HighOrderTetrahedronSetGeometryAlgorithms.inl"
#include <sofa/core/visual/VisualParams.h>
#include "BezierTetrahedronCoefficients.h"
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


double multinomial5Dim(const size_t n, sofa::defaulttype::Vec<5,int> tbiIn)
 {
	size_t i,ival;
	
	// divide n! with the largest of the multinomial coefficient
	std::sort(tbiIn.begin(),tbiIn.end());
	ival=1;
	for (i=n;i>tbiIn[4];--i){
		ival*=i;
	}
    return(((double)ival)/(sofa::helper::factorial(tbiIn[0])*sofa::helper::factorial(tbiIn[1])*sofa::helper::factorial(tbiIn[2])*sofa::helper::factorial(tbiIn[3])));
 }

double multinomial(const size_t n,const TetrahedronIndexVector tbiIn)
 {
	size_t i,ival;
	TetrahedronIndexVector tbi=tbiIn;
	// divide n! with the largest of the multinomial coefficient
	std::sort(tbi.begin(),tbi.end());
	ival=1;
	for (i=n;i>tbi[3];--i){
		ival*=i;
	}
    return(((double)ival)/(sofa::helper::factorial(tbi[0])*sofa::helper::factorial(tbi[1])*sofa::helper::factorial(tbi[2])));
 }
template< class DataTypes>
 BezierTetrahedronSetGeometryAlgorithms< DataTypes >::BezierTetrahedronSetGeometryAlgorithms() : 
HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>()
    {
    }
template< class DataTypes>
 void BezierTetrahedronSetGeometryAlgorithms< DataTypes >::init()
{
	HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::init();
	// recovers the pointer to a BezierTetrahedronSetTopologyContainer
	/// get the this->degree of the Bezier tetrahedron
	this->degree=this->container->getDegree();
	/// store the tetrahedron bezier index for each tetrahedron
	this->tbiArray=this->container->getTetrahedronIndexArray();
	/// compute the Bernstein coefficient for each control point in a tetrahedron
	bernsteinCoefficientArray.clear();
	bernsteinCoeffMap.clear();
	bernsteinCoefficientArray.resize(this->tbiArray.size());

	TetrahedronIndexVector tbi;
	/// precompute the factorial of the degree.
	for (size_t i=0;i<this->tbiArray.size();++i) {
		tbi=this->tbiArray[i];
		bernsteinCoefficientArray[i]=multinomial(this->degree,tbi); 
		bernsteinCoeffMap.insert(std::pair<TetrahedronIndexVector,Real>(tbi,(Real)bernsteinCoefficientArray[i]));
	}
	/// insert coefficient for the inferior degree
	HighOrderDegreeType i,j,k,l;
	for (i=0;i<=(this->degree-1);++i) {
		for (j=0;j<=(this->degree-i-1);++j) {
			for (k=0;k<=(this->degree-j-i-1);++k) {
				l=this->degree-1-i-j-k;
				tbi=TetrahedronIndexVector(i,j,k,l);
				bernsteinCoeffMap.insert(std::pair<TetrahedronIndexVector,Real>(tbi,(Real)multinomial(this->degree-1,tbi)));
			}
		}
	}

//	checkCoefficients();

}




template<class DataTypes>
typename DataTypes::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunction(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
{
	Real  val=pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
    typename std::map<TetrahedronIndexVector,Real>::iterator it=bernsteinCoeffMap.find(tbi);
	if (it!=bernsteinCoeffMap.end()) {
		val*=(*it).second;
		return(val);
	} else {
		val*=multinomial(tbi[0]+tbi[1]+tbi[2]+tbi[3],tbi);
		return(val);
	}
}
 
template<class DataTypes>
void BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeNodalValueDerivatives(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate,  const VecCoord& p, Deriv dpos[4])
{
	/// the 4 derivatives
	
	TetrahedronIndexVector tbi;
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	// initialize dpos
	for (j=0;j<4;++j) 
		dpos[j]=Deriv();
	bool isRational=this->container->isRationalSpline(tetrahedronIndex);
	if (isRational) {
		const HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Real weight=(Real)0.0f;
		Real dweight[4];
		dweight[0]=0.0f;dweight[1]=0.0f;dweight[2]=0.0f;dweight[3]=0.0f;

		Deriv pos;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			pos+=val*wa[indexArray[i]]*p[indexArray[i]];
			weight+=val*wa[indexArray[i]];
			Vec4 dval(0,0,0,0);
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
				dpos[j]+=dval[j]*wa[indexArray[i]]*p[indexArray[i]];
				dweight[j]+=dval[j]*wa[indexArray[i]];
			}
		}
		for (j=0;j<4;++j) {
			dpos[j]=(Real)this->degree*dpos[j]/weight-(Real)this->degree*(dweight[j]/(weight*weight))*pos;
		}
		
		
	} else {
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			Vec4 dval(0,0,0,0);
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
				dpos[j]+=dval[j]*p[indexArray[i]];
			}
		}

	}
}

template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real 
	 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate, const typename DataTypes::VecCoord& p)
 {
	/// the 3 derivatives
	Coord dpos[3];
	
	TetrahedronIndexVector tbi;
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	bool isRational=this->container->isRationalSpline(tetrahedronIndex);
	if (isRational) {
		const HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Real weight=(Real)0.0f;
        Real dweight[3]= {0,0,0};
		dweight[0]=dweight[1]=dweight[2]=(Real)0.0f;
		Coord pos;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			Vec4 dval(0,0,0,0);
			pos+=wa[indexArray[i]]*val*p[indexArray[i]];
			weight+=wa[indexArray[i]]*val;
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
			}
			for (j=0;j<3;++j) {
				dpos[j]+=(dval[j]-dval[3])*wa[indexArray[i]]*p[indexArray[i]];
				dweight[j]+=(dval[j]-dval[3])*wa[indexArray[i]];
			}
		}
		// computes the derivatives of the ratio of the 2 polynomial terms
		for (j=0;j<3;++j) {
			dpos[j]=dpos[j]/weight-(dweight[j]/(weight*weight))*pos;
		}
	}else {
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			Vec4 dval(0,0,0,0);
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
			}
			for (j=0;j<3;++j) {
				dpos[j]+=(dval[j]-dval[3])*p[indexArray[i]];
			}
		}
	}

	
	return(tripleProduct(dpos[0],dpos[1],dpos[2]));
 }

template<class DataTypes>
void BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeDeCasteljeauPoints(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate,  const VecCoord& p, Coord dpos[4])
{
	/// the 4 derivatives
	
	TetrahedronIndexVector tbi;
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	// initialize dpos
	for (j=0;j<4;++j) 
		dpos[j]=Coord();
	bool isRational=this->container->isRationalSpline(tetrahedronIndex);
	if (isRational) {
		const HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Real weight=(Real)0.0f;
		Real dweight[4];
		dweight[0]=0.0f;dweight[1]=0.0f;dweight[2]=0.0f;dweight[3]=0.0f;

		Coord pos;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			pos+=val*wa[indexArray[i]]*p[indexArray[i]];
			weight+=val*wa[indexArray[i]];
			Vec4 dval(0,0,0,0);
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
				dpos[j]+=dval[j]*wa[indexArray[i]]*p[indexArray[i]];
				dweight[j]+=dval[j]*wa[indexArray[i]];
			}
		}
		for (j=0;j<4;++j) {
			dpos[j]=dpos[j]/weight-(dweight[j]/(weight*weight))*pos;
		}
		
		
	} else {
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2])*pow(barycentricCoordinate[3],tbi[3]);
			Vec4 dval(0,0,0,0);
			for (j=0;j<4;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=3;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%4],tbi[(j+k)%4]));
				}
				dpos[j]+=dval[j]*p[indexArray[i]];
			}
		}

	}


}



 template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Vec4 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunctionDerivatives(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
 {
     Real  val=computeShapeFunction(tbi,barycentricCoordinate);
     Vec4 dval(0,0,0,0);
     for(unsigned i=0;i<4;++i)
         if(tbi[i] && barycentricCoordinate[i])
             dval[i]=(Real)tbi[i]*val/barycentricCoordinate[i];
     return dval;
 }

 template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeShapeFunctionHessian(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)
 {
     Vec4 dval = computeShapeFunctionDerivatives(tbi,barycentricCoordinate);
     Mat44 ddval;
     for(unsigned i=0;i<4;++i)
         if(barycentricCoordinate[i])
             for(unsigned j=0;j<4;++j)
             {
                 if(i==j) { if(tbi[i]>1) ddval[j][i]=((Real)tbi[i]-1.)*dval[j]/barycentricCoordinate[i]; }
                 else { if(tbi[i]) ddval[j][i]=(Real)tbi[i]*dval[j]/barycentricCoordinate[i]; }
             }
     return ddval;
 }
 //template<class DataTypes>
 //typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 //{
	// return(topology::binomialVector<4,typename DataTypes::Real>(tbi1,tbi2)/topology::binomial<typename DataTypes::Real>(2*this->degree,3));
 //}


//
//  template<class DataTypes>
// typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
// {
//	 Mat44 mat;
//	 size_t l,m,n;
//	 Real coeff=1.0/(binomial<typename DataTypes::Real>(2*this->degree-2,3));
//	 TetrahedronIndexVector tbi3,tbi4;
//
//	 for(l=0;l<4;++l) {
//		 if ((tbi1[l]>0) && (tbi2[l]>0)) {
//			 tbi3=tbi1;
//			 tbi3[l]-=1;
//			 tbi4=tbi2;
//			 tbi4[l]-=1;
//			 mat[l][l]=coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
//		 } else {
//			 mat[l][l]=(Real)0.0;
//		 }
//	 }
//	 for(l=0;l<6;++l) {
//		 m=edgesInTetrahedronArray[l][0];
//		 n=edgesInTetrahedronArray[l][1];
//		 if ((tbi1[m]>0) && (tbi2[n]>0)) {
//			 tbi3=tbi1;
//			 tbi3[m]-=1;
//			 tbi4=tbi2;
//			 tbi4[n]-=1;
//			 //	std::cerr<<" using direct with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
//			 mat[m][n]=coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
//		 } else
//			  mat[m][n]=0;
//		 if ((tbi1[n]>0) && (tbi2[m]>0)) {
//			 tbi3=tbi1;
//			 tbi3[n]-=1;
//			 tbi4=tbi2;
//			 tbi4[m]-=1;
//			 //		std::cerr<<" using transpose with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
//			 mat[n][m]==coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
//		 }
//	 }
//	 return(mat);
//}

   template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {
	
	 Real val=topology::binomialVector<4,typename DataTypes::Real>(tbi1,tbi2)/(topology::binomial<typename DataTypes::Real>(2*this->degree,3)*topology::binomial<typename DataTypes::Real>(this->degree,this->degree));
	 return(val);
 }
 template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
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
		 return(bezierTetrahedronDegree2AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 3:
		 return(bezierTetrahedronDegree3AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 4:
 		 return(bezierTetrahedronDegree4AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 5:
		 return(bezierTetrahedronDegree5AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 default:
		 {

			 Real val=computeAffineMassCoefficient(tbi1,tbi2);
			 return(val);
		 }
	 }
 }

 size_t multinomial5Vector(const TetrahedronIndexVector tbi[5]) 
{
	size_t i,j,n;
	size_t val=1;
	sofa::defaulttype::Vec<5,int> v;
	for (i=0;i<4;++i) {
		n=0;
		for(j=0;j<5;++j) {
			v[j]=tbi[j][i];
			n+=v[j];
		}
		val*=multinomial5Dim(n,v);
	}
	return(val);
}
 template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,
	 const size_t indr[3])
 {
	 // compute the mass coefficient in closed form using an horrible expression
	 Real val=0;
	 unsigned int permutation[6][4]= {{0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {2, 1, 0, 3}};
	 unsigned int triplets[8][3] = {{1, 2, 3}, {1, 2, 4}, {1, 4, 4}, {4, 2, 3}, {4, 4, 3}, {4, 4, 4}, {1, 4, 3}, {4, 2, 4}};
	 int tripletSign[8] ={1, -1, 1, -1, 1, -1, -1, 1}; 
	 int sign3[6] = {1, -1, -1, 1, 1, -1};

	 size_t i,j,k;
	 TetrahedronIndexVector tbir[3],tbiVec[5];
	 for (i=0;i<3;++i) {
		 tbir[i]=this->tbiArray[indr[i]];
	 }
	 tbiVec[3]=tbi1;
	 tbiVec[4]=tbi2;
	 unsigned int newIndex[3];
	 sofa::defaulttype::Vec<5,int> v;
	 v[0]=this->degree-1; v[1]=this->degree-1; v[2]=this->degree-1; v[3]=this->degree; v[4]=this->degree;
	 
	 Real factor=this->degree*this->degree*this->degree/(6*topology::binomial<typename DataTypes::Real>(5*this->degree-3,3)*multinomial5Dim(5*this->degree-3,v));

	 Real tmp;

	 for (i=0;i<6;++i) {

		 for(j=0;j<8;++j) {
			 for(k=0;k<3;++k) {
				 newIndex[k]=permutation[i][triplets[j][k]-1];
			 }
			 if (tbir[0][newIndex[0]]*tbir[1][newIndex[1]]*tbir[2][newIndex[2]]>0) {
				 for(k=0;k<3;++k) {
					 tbiVec[k]=tbir[k];
					 tbiVec[k][newIndex[k]]-=1;
				 }

				 val+=(Real)(tripletSign[j]*sign3[i]*(Real)multinomial5Vector(tbiVec));
			 }
		 }
	 }
	 val*=factor;
	 return(val);
 }

template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Real BezierTetrahedronSetGeometryAlgorithms<DataTypes>::getRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,
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
			 val=bezierTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
			 break;
		 case 3:
			 val=bezierTetrahedronDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
			 break;
		 }
	  } else {
		 val=computeRegularMassCoefficient(tbi1,tbi2,indr);
	  }
	return val;
 }


 template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::computeAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
 {
	  Mat44 mat;
	 size_t l,m,n;
	 Real coeff=this->degree*this->degree/(binomial<typename DataTypes::Real>(2*this->degree-2,3)*binomial<typename DataTypes::Real>(this->degree-1,this->degree-1));
	 TetrahedronIndexVector tbi3,tbi4;

	 for(l=0;l<4;++l) {
		 if ((tbi1[l]>0) && (tbi2[l]>0)) {
			 tbi3=tbi1;
			 tbi3[l]-=1;
			 tbi4=tbi2;
			 tbi4[l]-=1;
			 mat[l][l]=coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
		 } else {
			 mat[l][l]=(Real)0.0;
		 }
	 }
	 for(l=0;l<6;++l) {
		 m=edgesInTetrahedronArray[l][0];
		 n=edgesInTetrahedronArray[l][1];
		 if ((tbi1[m]>0) && (tbi2[n]>0)) {
			 tbi3=tbi1;
			 tbi3[m]-=1;
			 tbi4=tbi2;
			 tbi4[n]-=1;
			 //	std::cerr<<" using direct with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
			 mat[m][n]=coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
		 } else
			  mat[m][n]=0;
		 if ((tbi1[n]>0) && (tbi2[m]>0)) {
			 tbi3=tbi1;
			 tbi3[n]-=1;
			 tbi4=tbi2;
			 tbi4[m]-=1;
			 //		std::cerr<<" using transpose with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
			 mat[n][m]=coeff*binomialVector<4,typename DataTypes::Real>(tbi3,tbi4);
		 }
	 }
	 return(mat);
 }
  template<class DataTypes>
 typename BezierTetrahedronSetGeometryAlgorithms<DataTypes>::Mat44 BezierTetrahedronSetGeometryAlgorithms<DataTypes>::getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)
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
				 val[i][j]=bezierTetrahedronDegree2AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 3:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=bezierTetrahedronDegree3AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 4:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=bezierTetrahedronDegree4AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 5:
		 for(k=0,i=0;i<4;++i) {
			 for(j=0;j<4;++j,++k) {
				 val[i][j]=bezierTetrahedronDegree5AffineStiffnessCoefficient[index+k];
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
void  BezierTetrahedronSetGeometryAlgorithms<DataTypes>::checkCoefficients()
 { 
	 


	 // number of control points
	 size_t Nd=(this->degree+1)*(this->degree+2)*(this->degree+3)/6;
	 // number of unique entries in the symmetric mass matrix
	 size_t nbMassEntries=Nd*(Nd+1)/2;

	 if (this->degree<6) {
		 // check affine mass coefficients
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
					 assert(nbMassEntries==arraysize(bezierTetrahedronDegree2AffineMassCoefficient)); 
					 assert(fabs(bezierTetrahedronDegree2AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 3:
					 assert(nbMassEntries==arraysize(bezierTetrahedronDegree3AffineMassCoefficient)); 
					 assert(fabs(bezierTetrahedronDegree3AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 4:
					 assert(nbMassEntries==arraysize(bezierTetrahedronDegree4AffineMassCoefficient)); 
					 assert(fabs(bezierTetrahedronDegree4AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 5:
					 assert(nbMassEntries==arraysize(bezierTetrahedronDegree5AffineMassCoefficient)); 
					 assert(fabs(bezierTetrahedronDegree5AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 }

			 }
		 } 	 
	 }
	 if (this->degree<6) {
		 // check affine stiffness coefficients
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
					 assert(nbStiffnessEntries==arraysize(bezierTetrahedronDegree2AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=bezierTetrahedronDegree2AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 3:
					 assert(nbStiffnessEntries==arraysize(bezierTetrahedronDegree3AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=bezierTetrahedronDegree3AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 4:
					 assert(nbStiffnessEntries==arraysize(bezierTetrahedronDegree4AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=bezierTetrahedronDegree4AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 5:
					 assert(nbStiffnessEntries==arraysize(bezierTetrahedronDegree5AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<4;++k) {
						 for(l=0;l<4;++l,++m) {
							 val2[k][l]=bezierTetrahedronDegree5AffineStiffnessCoefficient[index+m];
						 }
					 }
					
					 assert((oneNorm(val-val2))<1e-8);
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
		 size_t aSize=arraysize(bezierTetrahedronDegree2RegularMassCoefficient);
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

							 Real val2=bezierTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)];
							 assert(fabs(bezierTetrahedronDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)*(Nd-r-2)/6+ (Nd-r-1)*(Nd-r-2)/2-(Nd-s-1)*(Nd-s)/2 +(t-s-1)]-val)<1e-8);
						 }
					 }
				 }
			 }
		 } 
	 }
	 //else {
		// 		 size_t i,j,r,s,t;
		// size_t Ndd=(this->degree+1)*(this->degree+2)*(this->degree)/6;
		// size_t nbElem=Nd*(Nd-1)*(Nd-2)/6;
		// size_t sss=this->tbiArray.size();
 	//	 std::vector<Real> coeffArray(nbElem*nbMassEntries);

		// for (i=0;i<this->tbiArray.size();++i) {

		//	 size_t ai=this->container->getHierarchicalIndex(i);
		//	 for (j=i;j<this->tbiArray.size();++j) {
		//		 size_t aj=this->container->getHierarchicalIndex(j);


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

		//  // write array in file
		// std::ostringstream st;
		// st<< "BezierTetrahedronDegree"<<(int)this->degree<<"RegMass.h";
		// std::string fname=st.str();

		// std::ofstream myfile("toto.h");
		// myfile.precision(16);
		// myfile <<"double bezierTetrahedronDegree"<<(int)this->degree<<"RegularMassCoefficient[]={";
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
