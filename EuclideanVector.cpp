//
// Created by Haoran Zhi on 16/8/16.
//

#include"EuclideanVector.h"
#include<iostream>
#include<cmath>
namespace evec {
    /* For all of my constructor, norm is initialised norm{-1} to indicate that norm has not been calculated ,
     * and every time that ecd_vector was changed , the value of norm will be set to -1 ,to indicate for next time
     * get the norm , the value of norm should be calculated again .
     * Otherwise the value of norm will be directly returned by the get norm
     */



    //default constructor with default value dim = 1, magnitude = 0.0
    EuclideanVector::EuclideanVector(unsigned int dim, double magnitude):
            dimensions{dim},
            ecd_vector{nullptr},
            norm{-1}
    {

        ecd_vector = new double[dimensions];

        for (unsigned int i = 0; i < dimensions; ++i) {
            ecd_vector[i] = magnitude;
        }
    }

    // Initializer_list constructor;
    EuclideanVector::EuclideanVector(const std::initializer_list<double>& lis):
            dimensions{1},
            ecd_vector{nullptr},
            norm{-1}
        {
            dimensions = lis.size();
            ecd_vector = new double[dimensions];

            unsigned int i =0;
            for(const auto &v:lis){
                ecd_vector[i++] = v;
            }
        }



    // copy constructor
    EuclideanVector::EuclideanVector(const EuclideanVector& rhs):
            dimensions{rhs.dimensions},
            ecd_vector{nullptr},
            norm{-1}
    {

        ecd_vector = new double[rhs.dimensions];
        for (unsigned int i = 0; i < rhs.dimensions; ++i) {
            ecd_vector[i] = rhs[i];
        }

    }


    // move constructor
    EuclideanVector::EuclideanVector(EuclideanVector&& rhs) noexcept :
            dimensions{std::move(rhs.dimensions)},
            ecd_vector{std::move(rhs.ecd_vector)},
            norm{std::move(rhs.norm)}
    {
        rhs.dimensions = 0;
        rhs.ecd_vector = nullptr;
        rhs.norm = -1;

    }

    // destructor
    EuclideanVector::~EuclideanVector() {
        //std::cout<<"destructor"<<'\n';
        delete[] ecd_vector;
    }

    //copy assignment
    EuclideanVector &EuclideanVector::operator=(const EuclideanVector& rhs) {
        // checking  self assignment
        if (*this != rhs) {
            delete[] ecd_vector;
            dimensions = rhs.dimensions;
            ecd_vector = new double[dimensions];
            for (unsigned int i = 0; i < dimensions; ++i) {
                ecd_vector[i] = rhs[i];
            }
            norm = rhs.norm;

        }
        return *this;

    }

    //move assignment
    EuclideanVector& EuclideanVector::operator=(EuclideanVector&& rhs) noexcept {
        //checking self assignment
        if (*this != rhs) {

            delete[] ecd_vector;


            dimensions = rhs.dimensions;
            ecd_vector = rhs.ecd_vector;
            norm = rhs.norm;


            // release the resources hold by rhs
            rhs.dimensions = 0;
            rhs.ecd_vector = nullptr;
            rhs.norm = -1;
        }
        return *this;
    }

    unsigned int EuclideanVector::getNumDimensions() const {
        return dimensions;
    }

    double EuclideanVector::get(unsigned int i) const {
        return ecd_vector[i];
    }

    double EuclideanVector::getEuclideanNorm() const {
        // norm == -1 indicates norm has not been calculated or
        // has been changed , therefore need to be calculated again.
        // and norm is mutable
        if (norm < 0) {
            double sum = 0.0;
            for (unsigned int i = 0; i < dimensions; ++i) {
                sum += pow(ecd_vector[i], 2);
            }
            norm = sqrt(sum);
        }

        return norm;

    }

    EuclideanVector EuclideanVector::createUnitVector() const {
        double div = getEuclideanNorm();
        EuclideanVector euc(dimensions, 0.0);
        // check whether div == 0
        // otherwise there will be problems
        if(div){
            for (unsigned int i = 0; i < dimensions; ++i) {
                euc[i] = ecd_vector[i] / div;
            }
        }
        return euc;
    }

    double EuclideanVector::operator[](unsigned int i) const {
        return ecd_vector[i];
    }

    // operator[] returns a reference to modify the value
    double& EuclideanVector::operator[](unsigned int i) {
        norm = -1;
        return ecd_vector[i];
    }


    EuclideanVector& EuclideanVector::operator+=(const EuclideanVector& rhs) {
        // only vectors with same dimensions will be added .
        if(dimensions == rhs.dimensions)
        {
            for (unsigned int i = 0; i < dimensions; ++i) {
                ecd_vector[i] += rhs[i];
            }

            norm = -1;
        }
        return *this;
    }

    EuclideanVector& EuclideanVector::operator-=(const EuclideanVector& rhs) {
        // only vectors with same dimensions will be operated .
        if(dimensions == rhs.dimensions)
        {
            for (unsigned int i = 0; i < dimensions; ++i) {
                ecd_vector[i] -= rhs[i];
            }

            norm = -1;
        }
        return *this;
    }

    EuclideanVector &EuclideanVector::operator*=(const int rhs) {
        for (unsigned int i = 0; i < dimensions; ++i) {
            ecd_vector[i] *= rhs;
        }

        norm = -1;
        return *this;
    }

    EuclideanVector &EuclideanVector::operator/=(const int rhs) {
        // check whether rhs == 0
        if (rhs) {
            for (unsigned int i = 0; i < dimensions; ++i) {
                ecd_vector[i] /= rhs;
            }

            norm = -1;
        }
        return *this;
    }

    // casting to vector
    EuclideanVector::operator std::vector<double>() const {
        std::vector<double> vec;
        for (unsigned int i = 0; i < dimensions; ++i) {
            vec.push_back(ecd_vector[i]);
        }

        return vec;
    }

    //casting to list
    EuclideanVector::operator std::list<double>() const {
        std::list<double> lis;
        for (unsigned int i = 0; i < dimensions; ++i) {
            lis.push_back(ecd_vector[i]);
        }
        return lis;
    }


    bool operator==(const EuclideanVector &lhs, const EuclideanVector &rhs) {

        /* Here are 3 return statements
         I do not know whether there are rules for the numbers of return statements
         I checked c++ core guideline , and the function return lists does not
         have specific rules for the number of this
         therefore could you please tell me whether there are problems with
         this
         */
        if ((lhs.dimensions) != (rhs.dimensions)) {
            return false;
        }
        else {
            for (unsigned int i = 0; i < lhs.dimensions; ++i) {
                if (lhs[i] != rhs[i]) {
                    return false;
                }
            }


            return true;
        }
    }

    bool operator!=(const EuclideanVector &lhs, const EuclideanVector &rhs) {
        return !(lhs == rhs);
    }


    EuclideanVector operator+(const EuclideanVector& lhs, const EuclideanVector& rhs) {

        EuclideanVector euc(lhs.dimensions, 0.0);
        // vectors with the dimensions be operated .
        // otherwise return an default vector
        if((lhs.dimensions) == (rhs.dimensions))
        {
            for (unsigned int i = 0; i < lhs.dimensions; ++i) {
                euc[i] = lhs[i] + rhs[i];
            }
        }
        return euc;
    }


    EuclideanVector operator-(const EuclideanVector &lhs, const EuclideanVector &rhs) {

        EuclideanVector euc(lhs.dimensions, 0.0);
        // vectors with the dimensions be operated .
        // otherwise return an default vector
        if((lhs.dimensions) == (rhs.dimensions))
        {
            for (unsigned int i = 0; i < lhs.dimensions; ++i) {
                euc[i] = lhs[i] - rhs[i];
            }
        }
        return euc;
    }

    double operator*(const EuclideanVector &lhs, const EuclideanVector &rhs) {
        double sum = 0.0;
        if((lhs.dimensions) == (rhs.dimensions))
        {
            for (unsigned int i = 0; i < lhs.dimensions; ++i) {
                sum += (lhs[i] * rhs[i]);
            }
        }
        return sum;
    }

    EuclideanVector operator*(const EuclideanVector &lhs, const int rhs) {

        EuclideanVector euc(lhs.dimensions, 0.0);
        for (unsigned int i = 0; i < lhs.dimensions; ++i) {
            euc[i] = lhs[i] * rhs;
        }
        return euc;
    }

    EuclideanVector operator*(const int lhs, const EuclideanVector &rhs) {

        EuclideanVector euc(rhs.dimensions, 0.0);
        for (unsigned int i = 0; i < rhs.dimensions; ++i) {
            euc[i] = rhs[i] * lhs;
        }
        return euc;
    }

    EuclideanVector operator/(const EuclideanVector &lhs, const int rhs) {
        // here check whether rhs ==0 , because 0 cannot be a divisior
        EuclideanVector euc(lhs.dimensions, 0.0);
        if (rhs) {

            for (unsigned int i = 0; i < lhs.dimensions; ++i) {
                euc[i] = lhs[i] / rhs;
            }


        }
        return euc;
    }

    // overloading the << operator
    // returns a reference to std::ostream
    // and if dimensions == 0
    // output a []
    std::ostream &operator<<(std::ostream& ofs, const EuclideanVector& rhs) {

        if (rhs.dimensions) {
            ofs << '[';
            unsigned int i = 0;
            for (i = 0; i < rhs.dimensions - 1; ++i) {
                ofs << rhs[i] << " ";
            }
            ofs << rhs[i] << ']';
        }
        else {
            ofs << '[' << ']';
        }
        return ofs;
    }
}


int main(){
evec::EuclideanVector a(3,5.0);
evec::EuclideanVector b(3,2.0);
	for(int i=0;i<10000;++i){
	a += b;
	a -=b ;
	a *= 5;
	a /= 5;
}
	std::cout<<a<<'\n';
	
}


