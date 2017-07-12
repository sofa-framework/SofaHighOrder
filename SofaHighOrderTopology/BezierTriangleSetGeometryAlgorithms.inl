#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_INL

#include "BezierTriangleSetGeometryAlgorithms.h"
#include "HighOrderTriangleSetGeometryAlgorithms.inl"
#include <sofa/core/visual/VisualParams.h>
#include "BezierTriangleCoefficients.h"
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/rmath.h>
namespace sofa
{

namespace component
{

namespace topology
{


const unsigned int edgesInTriangleArray[3][2] = {{0,1}, {0,2}, {1,2}};


double multinomial4Dim(const size_t n, sofa::defaulttype::Vec<4,int> tbiIn)
 {
	size_t i,ival;
	
	// divide n! with the largest of the multinomial coefficient
	std::sort(tbiIn.begin(),tbiIn.end());
	ival=1;
	for (i=n;i>tbiIn[3];--i){
		ival*=i;
	}
    return(((double)ival)/(sofa::helper::factorial(tbiIn[0])*sofa::helper::factorial(tbiIn[1])*sofa::helper::factorial(tbiIn[2])));
 }

double multinomial(const size_t n,const TriangleIndexVector tbiIn)
 {
	size_t i,ival;
	TriangleIndexVector tbi=tbiIn;
	// divide n! with the largest of the multinomial coefficient
	std::sort(tbi.begin(),tbi.end());
	ival=1;
	for (i=n;i>tbi[2];--i){
		ival*=i;
	}
    return(((double)ival)/(sofa::helper::factorial(tbi[0])*sofa::helper::factorial(tbi[1])));
 }
template< class DataTypes>
 BezierTriangleSetGeometryAlgorithms< DataTypes >::BezierTriangleSetGeometryAlgorithms() : 
HighOrderTriangleSetGeometryAlgorithms<DataTypes>()
    {
    }
template< class DataTypes>
 void BezierTriangleSetGeometryAlgorithms< DataTypes >::init()
{
	HighOrderTriangleSetGeometryAlgorithms<DataTypes>::init();
	// recovers the pointer to a BezierTriangleSetTopologyContainer
	/// get the this->degree of the Bezier triangle
	this->degree=this->container->getDegree();
	/// store the triangle bezier index for each triangle
	this->tbiArray=this->container->getTriangleIndexArray();
	/// compute the Bernstein coefficient for each control point in a triangle
	bernsteinCoefficientArray.clear();
	bernsteinCoeffMap.clear();
	bernsteinCoefficientArray.resize(this->tbiArray.size());

	TriangleIndexVector tbi;
	/// precompute the factorial of the degree.
	for (size_t i=0;i<this->tbiArray.size();++i) {
		tbi=this->tbiArray[i];
		bernsteinCoefficientArray[i]=multinomial(this->degree,tbi); 
		bernsteinCoeffMap.insert(std::pair<TriangleIndexVector,Real>(tbi,(Real) bernsteinCoefficientArray[i]));
	}
	/// insert coefficient for the inferior degree
	HighOrderDegreeType i,j,k,m;
	for (i=0;i<=(this->degree-1);++i) {
		for (j=0;j<=(this->degree-i-1);++j) {
			k=this->degree-1-i-j;
			tbi=TriangleIndexVector(i,j,k);
			bernsteinCoeffMap.insert(std::pair<TriangleIndexVector,Real>(tbi,(Real) multinomial(this->degree-1,tbi)));
		}
	}

//	checkCoefficients();

}




template<class DataTypes>
typename DataTypes::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunction(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
{
	Real  val=pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
    typename std::map<TriangleIndexVector,Real>::iterator it=bernsteinCoeffMap.find(tbi);
	if (it!=bernsteinCoeffMap.end()) {
		val*=(*it).second;
		return(val);
	} else {
        Real val2 = multinomial(this->degree, tbi);
        bernsteinCoeffMap.insert(std::pair<TriangleIndexVector, Real>(tbi, val2));
        return(val*val2);
	}
}
template<class DataTypes>
typename DataTypes::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionOfGivenDegree(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate, const HighOrderDegreeType deg)
{
    Real  val = pow(barycentricCoordinate[0], tbi[0])*pow(barycentricCoordinate[1], tbi[1])*pow(barycentricCoordinate[2], tbi[2]);
    typename std::map<TriangleIndexVector, Real>::iterator it = bernsteinCoeffMap.find(tbi);
    if (it != bernsteinCoeffMap.end()) {
        val *= (*it).second;
        return(val);
    }
    else {
        Real val2 = multinomial(deg , tbi);
        bernsteinCoeffMap.insert(std::pair<TriangleIndexVector, Real>(tbi, val2));
        return(val*val2);
    }
}
template<class DataTypes>
void BezierTriangleSetGeometryAlgorithms<DataTypes>::computeNodalValueDerivatives(const size_t triangleIndex, const Vec3 barycentricCoordinate,  const VecCoord& p, Deriv dpos[4])
{
	/// the 4 derivatives
	
	TriangleIndexVector tbi;
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	// initialize dpos
	for (j=0;j<3;++j) 
		dpos[j]=Deriv();
	bool isRational=this->container->isRationalSpline(triangleIndex);
	if (isRational) {
		const HighOrderTriangleSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Real weight=(Real)0.0f;
		Real dweight[3];
		dweight[0]=0.0f;dweight[1]=0.0f;dweight[2]=0.0f;

		Deriv pos;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			pos+=val*wa[indexArray[i]]*p[indexArray[i]];
			weight+=val*wa[indexArray[i]];
			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
				dpos[j]+=dval[j]*wa[indexArray[i]]*p[indexArray[i]];
				dweight[j]+=dval[j]*wa[indexArray[i]];
			}
		}
		for (j=0;j<3;++j) {
			dpos[j]=(Real)this->degree*dpos[j]/weight-(Real)this->degree*(dweight[j]/(weight*weight))*pos;
		}
		
		
	} else {
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];

				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
				dpos[j]+=dval[j]*p[indexArray[i]];
			}
		}

	}
}

template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Real 
	 BezierTriangleSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate, const typename DataTypes::VecCoord& p)
 {
	/// the 2 derivatives
	Deriv dpos[2];
	
	bool isRational=this->container->isRationalSpline(triangleIndex);
	TriangleIndexVector tbi;
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	if (isRational) {
		const HighOrderTriangleSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Real weight=(Real)0.0f;
        Real dweight[2]={0,0};
		Coord pos;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			// compute the numerator and denominator
			pos+=wa[indexArray[i]]*val*p[indexArray[i]];
			weight+=wa[indexArray[i]]*val;

			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
			}
			for (j=0;j<2;++j) {
				dpos[j]+=(dval[j]-dval[2])*wa[indexArray[i]]*p[indexArray[i]];
				dweight[j]+=(dval[j]-dval[2])*wa[indexArray[i]];
			}
		}
		// computes the derivatives of the ratio of the 2 polynomial terms
		for (j=0;j<2;++j) {
			dpos[j]=dpos[j]/weight-(dweight[j]/(weight*weight))*pos;
		}
	} else {
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
			}
			for (j=0;j<2;++j) {
				dpos[j]+=(dval[j]-dval[2])*p[indexArray[i]];
			}
		}
	}
	return(areaProduct(dpos[0],dpos[1]));
 }
template<class DataTypes>
void BezierTriangleSetGeometryAlgorithms<DataTypes>::computeDeCasteljeauPoints(const size_t triangleIndex, const Vec3 barycentricCoordinate,  const VecCoord& p, Coord dpos[3])
{
	
	/// the 3 derivatives

	TriangleIndexVector tbi;
	bool isRational=this->container->isRationalSpline(triangleIndex);
	size_t j;
	Real val;
	const VecPointID &indexArray=this->container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	// initialize dpos
	for (j=0;j<3;++j) 
		dpos[j]=Coord();
	if (isRational) {
		Real weight=(Real)0.0f;
		Real dweight[3];
		const HighOrderTriangleSetTopologyContainer::SeqWeights &wa=this->container->getWeightArray();
		Coord pos;
		dweight[0]=dweight[1]=0.0;
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			pos+=val*p[indexArray[i]];
			weight+=val*wa[indexArray[i]];
			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
				}else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
				dpos[j]+=dval[j]*p[indexArray[i]];
				dweight[j]+=dval[j]*wa[indexArray[i]];
			}
		}
		for (j=0;j<3;++j) {
			dpos[j]=dpos[j]/weight-(dweight[j]/(weight*weight))*pos;
		}
	} else{
		for(size_t i=0; i<this->tbiArray.size(); ++i)
		{
			tbi=this->tbiArray[i];
			val=bernsteinCoefficientArray[i]*pow(barycentricCoordinate[0],tbi[0])*pow(barycentricCoordinate[1],tbi[1])*pow(barycentricCoordinate[2],tbi[2]);
			Vec3 dval(0,0,0);
			for (j=0;j<3;++j) {
				if(tbi[j] && barycentricCoordinate[j]){
					dval[j]=(Real)tbi[j]*val/barycentricCoordinate[j];
					
				} else if ((barycentricCoordinate[j]==0.0f)&&(tbi[j]==1)) {
					dval[j]=bernsteinCoefficientArray[i];
					for (size_t k=1;k<=2;++k)
						dval[j]*=(Real)(pow( barycentricCoordinate[(j+k)%3],tbi[(j+k)%3]));
				}
				dpos[j]+=dval[j]*p[indexArray[i]];
			}
		}
	}

	
}



 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Vec3 BezierTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionDerivatives(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
 {
     Real  val=computeShapeFunction(tbi,barycentricCoordinate);
     Vec3 dval(0,0,0);
     for(unsigned i=0;i<3;++i)
         if(tbi[i] && barycentricCoordinate[i])
             dval[i]=(Real)tbi[i]*val/barycentricCoordinate[i];
     return dval;
 }
 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Vec3 BezierTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionDerivativesOfGivenDegree(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate, const HighOrderDegreeType deg)
 {
     Real  val = computeShapeFunctionOfGivenDegree(tbi, barycentricCoordinate, deg);
     Vec3 dval(0, 0, 0);
     for (unsigned i = 0; i<3; ++i)
         if (tbi[i] && barycentricCoordinate[i])
             dval[i] = (Real)tbi[i] * val / barycentricCoordinate[i];
     return dval;
 }
 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Mat33 BezierTriangleSetGeometryAlgorithms<DataTypes>::computeShapeFunctionHessian(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)
 {
     Vec3 dval = computeShapeFunctionDerivatives(tbi,barycentricCoordinate);
     Mat33 ddval;
     for(unsigned i=0;i<3;++i)
         if(barycentricCoordinate[i])
             for(unsigned j=0;j<3;++j)
             {
                 if(i==j) { if(tbi[i]>1) ddval[j][i]=((Real)tbi[i]-1.)*dval[j]/barycentricCoordinate[i]; }
                 else { if(tbi[i]) ddval[j][i]=(Real)tbi[i]*dval[j]/barycentricCoordinate[i]; }
             }
     return ddval;
 }


   template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::computeAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	
	 Real val=topology::binomialVector<3,typename DataTypes::Real>(tbi1,tbi2)/(topology::binomial<typename DataTypes::Real>(2*this->degree,2)*topology::binomial<typename DataTypes::Real>(this->degree,this->degree));
	 return(val);
 }
 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::getAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	 size_t i=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi1[0])*(this->degree+2-tbi1[0])/2+  tbi1[1]; 
	 size_t j=  (this->degree+1)*(this->degree+2)/2-(this->degree+1-tbi2[0])*(this->degree+2-tbi2[0])/2+  tbi2[1]; 
	 size_t ii=std::min(i,j);
     size_t jj=std::max(i,j);
	 size_t s=this->tbiArray.size();
	 switch(this->degree) 
	 {
	 case 2:
		 return(bezierTriangleDegree2AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 3:
		 return(bezierTriangleDegree3AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 4:
 		 return(bezierTriangleDegree4AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 case 5:
		 return(bezierTriangleDegree5AffineMassCoefficient[ ii*(2*s-ii+1)/2+(jj-ii)] );
		 break;
	 default:
		 {

			 Real val=computeAffineMassCoefficient(tbi1,tbi2);
			 return(val);
		 }
	 }
 }

 size_t multinomial4Vector(const TriangleIndexVector tbi[4]) 
{
	size_t i,j,n;
	size_t val=1;
	sofa::defaulttype::Vec<4,int> v;
	for (i=0;i<3;++i) {
		n=0;
		for(j=0;j<4;++j) {
			v[j]=tbi[j][i];
			n+=v[j];
		}
		val*=multinomial4Dim(n,v);
	}
	return(val);
}
 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::computeRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,
	 const size_t indr[2])
 {
	 // compute the mass coefficient in closed form using an horrible expression
	 Real val=0;
	 unsigned int permutation[2][3]= {{0, 1, 2}, {1, 0, 2}};
	 unsigned int triplets[4][2] = {{1, 2}, {1, 3}, {3, 2}, {3, 3}};
	 int tripletSign[4] ={1, -1, -1, 1}; 
	 int sign3[2] = {1, -1};

	 size_t i,j,k;
	 TriangleIndexVector tbir[2],tbiVec[4];
	 for (i=0;i<2;++i) {
		 tbir[i]=this->tbiArray[indr[i]];
	 }
	 tbiVec[2]=tbi1;
	 tbiVec[3]=tbi2;
	 unsigned int newIndex[2];
	 sofa::defaulttype::Vec<4,int> v;
	 v[0]=this->degree-1; v[1]=this->degree-1; v[2]=this->degree-1; v[3]=this->degree; 
	 
	 Real factor=this->degree*this->degree*this->degree/(2*topology::binomial<typename DataTypes::Real>(4*this->degree-2,2)*multinomial4Dim(4*this->degree-2,v));


	 for (i=0;i<2;++i) {

		 for(j=0;j<4;++j) {
			 for(k=0;k<2;++k) {
				 newIndex[k]=permutation[i][triplets[j][k]-1];
			 }
			 if (tbir[0][newIndex[0]]*tbir[1][newIndex[1]]>0) {
				 for(k=0;k<2;++k) {
					 tbiVec[k]=tbir[k];
					 tbiVec[k][newIndex[k]]-=1;
				 }

				 val+=(Real)(tripletSign[j]*sign3[i]*(Real)multinomial4Vector(tbiVec));
			 }
		 }
	 }
	 val*=factor;
	 return(val);
 }

template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Real BezierTriangleSetGeometryAlgorithms<DataTypes>::getRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,
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
			 val=bezierTriangleDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];

			 break;
		 case 3:
			 val=bezierTriangleDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 case 4:
			 val=bezierTriangleDegree4RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 case 5:
			 val=bezierTriangleDegree5RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
			 break;
		 }


	 } else {
		 val=computeRegularMassCoefficient(tbi1,tbi2,indr);
	 }
	 return val;
 }


 template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Mat33 BezierTriangleSetGeometryAlgorithms<DataTypes>::computeAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
 {
	  Mat33 mat;
	 size_t l,m,n;
	 Real coeff=this->degree*this->degree/(binomial<typename DataTypes::Real>(2*this->degree-2,2)*binomial<typename DataTypes::Real>(this->degree-1,this->degree-1));
	 TriangleIndexVector tbi3,tbi4;

	 for(l=0;l<3;++l) {
		 if ((tbi1[l]>0) && (tbi2[l]>0)) {
			 tbi3=tbi1;
			 tbi3[l]-=1;
			 tbi4=tbi2;
			 tbi4[l]-=1;
			 mat[l][l]=coeff*binomialVector<3,typename DataTypes::Real>(tbi3,tbi4);
		 } else {
			 mat[l][l]=(Real)0.0;
		 }
	 }
	 for(l=0;l<3;++l) {
		 m=edgesInTriangleArray[l][0];
		 n=edgesInTriangleArray[l][1];
		 if ((tbi1[m]>0) && (tbi2[n]>0)) {
			 tbi3=tbi1;
			 tbi3[m]-=1;
			 tbi4=tbi2;
			 tbi4[n]-=1;
			 //	std::cerr<<" using direct with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
			 mat[m][n]=coeff*binomialVector<3,typename DataTypes::Real>(tbi3,tbi4);
		 } else
			  mat[m][n]=0;
		 if ((tbi1[n]>0) && (tbi2[m]>0)) {
			 tbi3=tbi1;
			 tbi3[n]-=1;
			 tbi4=tbi2;
			 tbi4[m]-=1;
			 //		std::cerr<<" using transpose with tbi1="<< Vec<4,unsigned int>(tbi1)<< "tbi2="<<Vec<4,unsigned int>(tbi2)<<"tbi3="<<Vec<4,unsigned int>(tbi3)<<" and tbi4="<<Vec<4,unsigned int>(tbi4)<<std::endl;
			 mat[n][m]=coeff*binomialVector<3,typename DataTypes::Real>(tbi3,tbi4);
		 }
	 }
	 return(mat);
 }
  template<class DataTypes>
 typename BezierTriangleSetGeometryAlgorithms<DataTypes>::Mat33 BezierTriangleSetGeometryAlgorithms<DataTypes>::getAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)
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
				 val[i][j]=bezierTriangleDegree2AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 3:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=bezierTriangleDegree3AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 4:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=bezierTriangleDegree4AffineStiffnessCoefficient[index+k];
			 }
		 }
		 break;
	 case 5:
		 for(k=0,i=0;i<3;++i) {
			 for(j=0;j<3;++j,++k) {
				 val[i][j]=bezierTriangleDegree5AffineStiffnessCoefficient[index+k];
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
void  BezierTriangleSetGeometryAlgorithms<DataTypes>::checkCoefficients()
 { 
	 


	 // number of control points
	 size_t Nd=(this->degree+1)*(this->degree+2)/2;
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
				 TriangleIndexVector tbi1,tbi2;
				 tbi1=this->tbiArray[aai];
				 tbi2=this->tbiArray[aaj];
				 Real val=computeAffineMassCoefficient(tbi1,tbi2);

				 switch(this->degree) {
				 case 2:
					 assert(nbMassEntries==arraysize(bezierTriangleDegree2AffineMassCoefficient)); 
					 assert(fabs(bezierTriangleDegree2AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 3:
					 assert(nbMassEntries==arraysize(bezierTriangleDegree3AffineMassCoefficient)); 
					 assert(fabs(bezierTriangleDegree3AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 4:
					 assert(nbMassEntries==arraysize(bezierTriangleDegree4AffineMassCoefficient)); 
					 assert(fabs(bezierTriangleDegree4AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 case 5:
					 assert(nbMassEntries==arraysize(bezierTriangleDegree5AffineMassCoefficient)); 
					 assert(fabs(bezierTriangleDegree5AffineMassCoefficient[ i*(2*s-i+1)/2+(j-i)]-val)<1e-8);
					 break;
				 }

			 }
		 } 	 
	 }
	 if (this->degree<6) {
		 // check affine stiffness coefficients
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
					 assert(nbStiffnessEntries==arraysize(bezierTriangleDegree2AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=bezierTriangleDegree2AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 3:
					 assert(nbStiffnessEntries==arraysize(bezierTriangleDegree3AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=bezierTriangleDegree3AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 4:
					 assert(nbStiffnessEntries==arraysize(bezierTriangleDegree4AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=bezierTriangleDegree4AffineStiffnessCoefficient[index+m];
						 }
					 }
					 assert((oneNorm(val-val2))<1e-8);
					 break;
				 case 5:
					 assert(nbStiffnessEntries==arraysize(bezierTriangleDegree5AffineStiffnessCoefficient)); 
					 for(m=0,k=0;k<3;++k) {
						 for(l=0;l<3;++l,++m) {
							 val2[k][l]=bezierTriangleDegree5AffineStiffnessCoefficient[index+m];
						 }
					 }
					
					 assert((oneNorm(val-val2))<1e-8);
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
			  aSize=arraysize(bezierTriangleDegree2RegularMassCoefficient);
			  break;
		 case 3:
			  aSize=arraysize(bezierTriangleDegree3RegularMassCoefficient);
			  break;
		 case 4:
			  aSize=arraysize(bezierTriangleDegree4RegularMassCoefficient);
			  break;
		 case 5:
			  aSize=arraysize(bezierTriangleDegree5RegularMassCoefficient);
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
							 val2=bezierTriangleDegree2RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];

							 break;
						 case 3:
							 val2=bezierTriangleDegree3RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 case 4:
							 val2=bezierTriangleDegree4RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 case 5:
							 val2=bezierTriangleDegree5RegularMassCoefficient[ massIndex+nbElem-(Nd-r)*(Nd-r-1)/2+ (s-r-1)];
							 break;
						 }
						 assert(fabs(val2-val)<1e-8);
					 }
				 }

			 }
		 } 
	 }
}


} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_TETEAHEDRONSETGEOMETRYALGORITHMS_INL
