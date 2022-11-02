#include <iostream>
#include <cmath>
#include <cassert>
#pragma once

namespace math
{
    template<typename T>
    T* zeros(size_t N)
    {
        T* d = new T[N];
        for (size_t i = 0; i < N; ++i) d[i] = (T)(0);
        return d;
    }

    template<typename T>
    class matrix
    {
    protected:
        size_t N, M;
        T* data;

    public:
        // Constructor Routine
        matrix(): data(nullptr), N(0), M(0) {}
        matrix(const T& d): data(new T[1]), N(1), M(1) {data[0]=d;}
        matrix(size_t N, size_t M): data(new T[N*M]), N(N), M(M) {}
        matrix(size_t N): data(new T[N*N]), N(N), M(N){}
        matrix(T* data_, size_t N, size_t M): data(data_), N(N), M(M) {}
        size_t height() const {return N;}
        size_t width () const {return M;}
        size_t size  () const {return N*M;}

        //RAII
        matrix(const matrix& lhs): N(lhs.N), M(lhs.M)
        {
            data = new T[N*M];
            for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < M; ++j)
                data[i * M + j] = lhs.data[i * M + j];
        }
        matrix& operator=(const matrix& lhs)
        {
            if (this == &lhs)
                return *this;
            matrix<T> t(lhs);
            std::swap(data, t.data);
            std::swap(N, t.N);
            std::swap(M, t.M);
            return *this;
        }
        matrix(matrix&& rhs): data(rhs.data), N(N), M(M)
        {
            rhs.data = nullptr;
            rhs.N = 0;
            rhs.M = 0;
        }
        matrix& operator=(matrix&& rhs)
        {
            matrix<T> t(std::move(rhs));
            std::swap(data, t.data);
            std::swap(N, t.N);
            std::swap(M, t.M);

            return *this;
        }
        ~matrix()
        {
            delete[] data;
        }
        T* operator[](size_t i)
        {
            return data + (i * M);
        }
        T operator()(size_t i, size_t j) const
        {
            assert(i < N && j < M && i >= 0 && j >= 0);
            return data[i * M + j];
        }
        T& operator()(size_t i, size_t j)
        {
            assert(i < N && j < M && i >= 0 && j >= 0);
            return data[i * M + j];
        }
        matrix& operator=(T const &t)
        {
            for(size_t i = 0; i < N * M; ++i)
                data[i] = t;
            return *this;
        }
        
        // Operators Routine
        /*
            Сделать все бинарные операторы как +
        */
        bool    operator ==(const matrix other)
        {
            if (N != other.N || M != other.M) return false;
            bool f = true;
            for(size_t i = 0; i < N; ++i) for(size_t j = 0; j < M; ++j)
                f &= operator()(i,j) == other.operator()(i,j);
            return f;
        }
        bool    operator !=(const matrix other)
        {
            return !operator==(other);
        }
        matrix& operator  +(const matrix& other)
        {
            assert(N == other.N && M == other.M);
            matrix* summ = new matrix(N,M);
            for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < M; ++j)
                summ->operator()(i,j) = this->operator()(i,j) + other.operator()(i,j);
            return *summ;
        }
        matrix& operator  -(const matrix& other)
        {
            assert(N == other.N && M == other.M);
            matrix* summ = new matrix(N,M);
            for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < M; ++j)
                summ->operator()(i,j) = this->operator()(i,j) - other.operator()(i,j);
            return *summ;
        }
        matrix& operator  *(const matrix& other)
        {
            assert(M == other.N);
            matrix* mult = new matrix(N, other.M);
            *mult = (T)(0);
            for(size_t i = 0; i < mult->N; ++i) for(size_t j = 0; j < mult->M; ++j) for(size_t k = 0; k < M; ++k)
                mult->operator()(i,j) += operator()(i,k)*other.operator()(k,j);
            return *mult;
        }
        matrix& operator  /(const T d)
        {
            assert(d != (T)(0));
            matrix* m = new matrix(N,M);
            for(size_t i = 0; i < m->N; ++i) for(size_t j = 0; j < m->M; ++j)
                m->operator()(i,j) = operator()(i,j)/d;
            return *m;
        }
        matrix& operator +=(const matrix& other)
        {
            assert(N == other.N && M == other.M);
            for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < M; ++j)
                this->operator()(i,j) += other.operator()(i,j);
            return *this;
        }
        matrix& operator -=(const matrix& other)
        {
            assert(N == other.N && M == other.M);
            for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < M; ++j)
                this->operator()(i,j) -= other.operator()(i,j);
            return *this;
        }
        matrix& operator *=(const matrix& other)
        {
            assert(M == other.N);
            matrix mult(N, other.M); mult = (T)(0);
            for(size_t i = 0; i < mult.N; ++i)
                for(size_t j = 0; j < mult.M; ++j)
                    for(size_t k = 0; k < M; ++k)
                        mult.operator()(i,j) += this->operator()(i,k)*other.operator()(k,j);
            *this = mult;
            return *this;
        }
        matrix& operator /=(const T d)
        {
            assert(d != (T)(0));
            matrix<T>m(M,M); m = (T)(0);
            for (size_t i = 0; i < m.N; ++i) m(i,i) = 1/d;
            return this->operator*=(m);
        }
        matrix& transpose()
        {
            matrix m(M,N);
            for(size_t i = 0; i < M; ++i) for(size_t j = 0; j < N; ++j)
                m(i,j) = operator()(j,i);
            return m;
        }
        matrix& column(const size_t j)
        {
            assert(j >= 0 && j < M);
            matrix* sub = new matrix(N,1);
            for(size_t i = 0; i < N; ++i)
                sub->operator()(i,0) = operator()(i,j);
            return *sub;
        }
        matrix& row   (const size_t i)
        {
            assert(i >= 0 && i < N);
            matrix* sub = new matrix(1,M);
            for(size_t j = 0; j < M; ++j) 
                sub->operator()(0,j) = operator()(i,j);
            return *sub;
        }
        T conv()
        {
            assert(N == 1 && M == 1);
            return operator()(0,0);
        }
        void print()
        {
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < M; ++j) std::cout << operator()(i,j) << ' ';
                std::cout << std::endl;
            }
            return;
        }
    };

    template<typename T>
    matrix<T>& diag(T d, size_t N)
    {
        matrix<T>* m = new matrix<T>(N,N);
        *m = 0;
        for (size_t i = 0; i < m->height(); ++i)
            (*m)(i,i) = d;
        return *m;
    }
    template<typename T>
    matrix<T>& Id(size_t N)
    {
        return diag<T>((T)(1), N);
    }
    template<typename T>
    matrix<T>& One(size_t i, size_t j, size_t N)
    {
        assert(i >= 0 && j >= 0 && i <= N && j <= N);
        matrix<T>* m = new matrix(N,N); 
        (*m)(i, j) = (T)(1);
        return *m;
    }


    template<typename T>
    class spatial_vector: public matrix<T>
    {
    public:
        T x, y, z;

        private:
        void correct_xyz_parameters()
        {
            x = this->operator()(0);
            y = this->operator()(1);
            z = this->operator()(2);
            return;
        }

    public:
        spatial_vector(): matrix<T>(3,1)
        {
            this->operator()(0) = 0; x = 0;
            this->operator()(1) = 0; y = 0;
            this->operator()(2) = 0; z = 0;
        }
        spatial_vector(T x_, T y_, T z_): matrix<T>(3,1)
        {
            this->operator()(0) = x_; x = x_;
            this->operator()(1) = y_; y = y_;
            this->operator()(2) = z_; z = z_;
        }
        spatial_vector(matrix<T> m): matrix<T>(3,1)
        {
            assert(m.size() == 3);
            if (m.N == 1)
            {
                this->operator()(0,0) = m.operator()(0,0); x = m.operator()(0,0);
                this->operator()(0,1) = m.operator()(0,1); y = m.operator()(0,1);
                this->operator()(0,2) = m.operator()(0,2); x = m.operator()(0,2);
            }
            if (m.M == 1)
            {
                this->operator()(0,0) = m.operator()(0,0); x = m.operator()(0,0);
                this->operator()(1,0) = m.operator()(1,0); y = m.operator()(1,0);
                this->operator()(2,0) = m.operator()(2,0); x = m.operator()(2,0);
            }
            
        }
        T  operator()(size_t i) const
        {
            assert(i >= 0 && i < 3);
            return this->operator()(i,0);
        }
        T& operator()(size_t i)
        {
            assert(i >= 0 && i < 3);
            return this->data[i];
        }
        T  operator[](size_t i) const
        {
            assert(i >= 0 && i < 3);
            return this->data[i];
        }
        T& operator[](size_t i)
        {
            assert(i >= 0 && i < 3);
            return this->data[i];
        }
        void operator +=(spatial_vector other)
        {
            this->operator()(0) += other.operator()(0); 
            this->operator()(1) += other.operator()(1);
            this->operator()(2) += other.operator()(2);
            this->correct_xyz_parameters();
            return;
        }
        void operator -=(spatial_vector other)
        {
            this->operator()(0) -= other.operator()(0); 
            this->operator()(1) -= other.operator()(1);
            this->operator()(2) -= other.operator()(2);
            correct_xyz_parameters();
            return;
        }
        void operator *=(T d)
        {
            this->operator()(0) *= d;
            this->operator()(1) *= d;
            this->operator()(2) *= d;
            correct_xyz_parameters();
        }
        void operator /=(T d)
        {
            this->operator()(0) /= d;
            this->operator()(1) /= d;
            this->operator()(2) /= d;
            correct_xyz_parameters();
        }
        T    operator  *(spatial_vector other)
        {
            return ((matrix(this->data,3,1).transpose()).operator*(other)).conv();
        }
        T norm()
        {
            return sqrt(operator*(*this));
        }
        spatial_vector operator^(spatial_vector other)
        {
            return spatial_vector
            (
                this->operator()(1)*other.operator()(2) - this->operator()(2)*other.operator()(1),
                this->operator()(2)*other.operator()(0) - this->operator()(0)*other.operator()(2),
                this->operator()(0)*other.operator()(1) - this->operator()(1)*other.operator()(0)
            );
        }
        void operator^=(spatial_vector other)
        {
            this->operator()(0) = this->operator()(1)*other.operator()(2) - this->operator()(2)*other.operator()(1);
            this->operator()(1) = this->operator()(2)*other.operator()(0) - this->operator()(0)*other.operator()(2);
            this->operator()(2) = this->operator()(0)*other.operator()(1) - this->operator()(1)*other.operator()(0);
            correct_xyz_parameters();
        }
        void print()
        {
            std::cout << '(' << this->operator()(0) << ", " << 
                                this->operator()(1) << ", " <<
                                this->operator()(2) << ')' << std::endl;
        }
    };
    /*
        сделать так, чтобы оператор + у векторов возвращал вектора
    *//*
    class quaternion
    {
        public:

        double s;
        spatial_vector<double> v;

        quaternion(double a, double b, double c, double d): s(a), v(b, c, d) {}
        quaternion(double a, spatial_vector<double> v_): s(a), v(v_) {}
        quaternion(double a, matrix<double> m): s(a), v(m) {}

        void operator  =(quaternion other)
        {
            s = other.s; v = other.v;
            return;
        }
        bool operator ==(quaternion other)
        {
            if (s != other.s) return false;
            if (v != other.v) return false;
            return true;
        }
        bool operator !=(quaternion other)
        {
            return !(operator==(other));
        }
        quaternion operator  +(quaternion other)
        {
            //std::cout << s+other.s << std::endl; (v+other.v).print();
            //(quaternion(s+other.s, v+other.v)).print();
            (v+other.v).print();
            return quaternion(s+other.s, v+other.v);
        }
        quaternion operator  -(quaternion other)
        {
            return quaternion(s-other.s, v-other.v);
        }
    
        void print()
        {
            std::cout << 
                s << " + " << 
                v.x << " i + " <<
                v.y << " j + " <<
                v.z << " k" << std::endl;
            return;
        }
    };*/

};