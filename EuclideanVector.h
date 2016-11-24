//
// Created by Haoran Zhi on 16/8/16.
//

#ifndef INC_6771ASSIGNMEN2_EUCLIDEANVECTOR_H
#define INC_6771ASSIGNMEN2_EUCLIDEANVECTOR_H
#include<vector>
#include<list>
#include<functional>
#include <iterator>
#include<iostream>
#include<initializer_list>
namespace evec{
    class EuclideanVector {
        friend bool operator ==(const EuclideanVector&, const EuclideanVector&);
        friend bool operator !=(const EuclideanVector&, const EuclideanVector&);
        friend EuclideanVector operator+(const EuclideanVector&, const EuclideanVector&);
        friend EuclideanVector operator-(const EuclideanVector&, const EuclideanVector&);
        friend double operator*(const EuclideanVector&, const EuclideanVector&);
        friend EuclideanVector operator*(const EuclideanVector&, const int);
        friend EuclideanVector operator*(const int, const EuclideanVector&);
        friend EuclideanVector operator/(const EuclideanVector&, const int);
        friend std::ostream& operator<<(std::ostream&,const EuclideanVector&);
    public:

        EuclideanVector(unsigned int dim = 1,double magnitude = 0.0);

        //This is a template definition for using iterator as parameters.
        // the std::enable_if<> means if the expression in <> is true, then the template is used
        // otherwise the template is ignored .

        template<class InputIterator>
        EuclideanVector(const InputIterator& begin,const InputIterator& end,
                        typename std::enable_if<!std::is_integral<InputIterator>::value>::type* = 0):
                dimensions{0},
                ecd_vector{nullptr},
                norm{-1}
        {

            for(auto it = begin; it !=end; ++ it){

                ++dimensions;
            }


            ecd_vector = new double[dimensions];
            unsigned int i = 0;
            // Here i++ is used instead of ++i
            // is because the value of i(before adding) is useful
            for(auto it = begin; it!= end;++it){
                ecd_vector[i++] = *it;
            }
        };


        EuclideanVector( const std::initializer_list<double>& lis);

        // copy consturctor
        EuclideanVector(const EuclideanVector &);
        // move constructor
        EuclideanVector(EuclideanVector &&) noexcept ;

        ~EuclideanVector();
        // copy assignment
        EuclideanVector &operator=(const EuclideanVector&);
        // move assignment
        EuclideanVector &operator=(EuclideanVector &&)noexcept ;


        unsigned int getNumDimensions() const ;
        double get(unsigned int) const;
        double getEuclideanNorm() const;



        EuclideanVector createUnitVector() const;
        double operator[](unsigned int) const;
        double& operator[](unsigned int);


        EuclideanVector& operator+=(const EuclideanVector&);
        EuclideanVector& operator-=(const EuclideanVector&);
        EuclideanVector& operator*=(const int);
        EuclideanVector& operator/=(const int);


        operator std::vector<double>() const;
        operator std::list<double>() const;





    private:
        unsigned int dimensions;
        // the allocated memory to store double values
        double* ecd_vector;
        mutable double norm;
    };


    bool operator ==(const EuclideanVector&, const EuclideanVector&);
    bool operator !=(const EuclideanVector&, const EuclideanVector&);
    EuclideanVector operator+(const EuclideanVector&, const EuclideanVector&);
    EuclideanVector operator-(const EuclideanVector&, const EuclideanVector&);
    double operator*(const EuclideanVector&, const EuclideanVector&);
    EuclideanVector operator*(const EuclideanVector&, const int);
    EuclideanVector operator*(const int, const EuclideanVector&);
    EuclideanVector operator/(const EuclideanVector&, const int);
    std::ostream& operator<<(std::ostream& ,const EuclideanVector&);


}

#endif //INC_6771ASSIGNMEN2_EUCLIDEANVECTOR_H
