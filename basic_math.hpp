#include <iostream>
#include <cmath>
#include <cassert>

namespace math
{
    template<typename T>
    T* zeros(std::size_t N)
    {
        T* d = new T[N];
        for (std::size_t i = 0; i < N; ++i) d[i] = (T)(0);
        return d;
    }

    template<typename T>
    class matrix
    {
        protected:
        std::size_t N, M;
        T* data;

        public:

        matrix(T d): N(1), M(1) {data[0]=d;}
        matrix(std::size_t N, std::size_t M): data(zeros<T>(N*M)), N(N), M(M){}
        matrix(std::size_t N): data(zeros<T>(N*N)), N(N), M(N){}
        matrix(T* data_, std::size_t N, std::size_t M): data(data_), N(N), M(M) {}
        std::size_t height() {return this->N;}
        std::size_t width () {return this->M;}
        std::size_t size  () {return this->N*this->M;}
        T&     operator ()(std::size_t i, std::size_t j)
        {
            assert(i < this->N && j < this->M && i >= 0 && j >= 0);
            return this->data[i*this->M+j];
        }
        T&     operator ()(std::size_t i)
        {
            assert(this->N == 1 || this->M == 1);
            if (this->N == 1)
            {
                assert(i >= 0 && i < this->M);
                return this->operator()(0,i);
            }
            if (this->M == 1)
            {
                assert(i >= 0 && i < this->N);
                return this->operator()(i,0);
            }
            return this->operator()(0,0);
        }
        bool   operator ==(matrix other)
        {
            if (this->N != other.N || this->M != other.M) return false;
            bool f = true;
            for(std::size_t i = 0; i < this->N; ++i) for(std::size_t j = 0; j < this->M; ++j)
                f &= this->operator()(i,j) == other.operator()(i,j);
            return f;
        }
        bool   operator !=(matrix other)
        {
            return !this->operator==(other);
        }
        void   operator  =(matrix other)
        {
            this->N=other.N; this->M=other.M;
            for (std::size_t i = 0; i < this->N*this->M; ++i) this->data[i]=other.data[i];
            return;
        }
        matrix operator  +(matrix other)
        {
            assert(this->N == other.N && this->M == other.M);
            matrix summ(this->N, this->M);
            for (std::size_t i = 0; i < this->N; ++i) for (std::size_t j = 0; j < this->M; ++j)
                summ.operator()(i,j)=this->operator()(i,j) + other.operator()(i,j);
            return summ;
        }
        matrix operator  -(matrix other)
        {
            assert(this->N == other.N && this->M == other.M);
            matrix summ(this->N, this->M);
            for (std::size_t i = 0; i < this->N; ++i) for (std::size_t j = 0; j < this->M; ++j)
                summ.operator()(i,j)=this->operator()(i,j) - other.operator()(i,j);
            return summ;
        }
        matrix operator  *(matrix other)
        {
            assert(this->M == other.N || (this->size() == 1) || (other.size() == 1));
            if (this->size() == 1)
            {
                matrix m(other.N, other.M);
                for(std::size_t i = 0; i < m.N; ++i) for(std::size_t j = 0; j < m.M; ++j) 
                    m.operator()(i,j) = this->operator()(0,0)*other.operator()(i,j);
                return m;
            }
            if (other.size() == 1)
            {
                matrix m(this->N, this->M);
                for(std::size_t i = 0; i < m.N; ++i) for(std::size_t j = 0; j < m.M; ++j) 
                    m.operator()(i,j) = other.operator()(0,0)*this->operator()(i,j);
                return m;
            }
            matrix mult(this->N, other.M);
            for(std::size_t i = 0; i < mult.N; ++i) for(std::size_t j = 0; j < mult.M; ++j) for(std::size_t k = 0; k < this->M; ++k)
                mult.operator()(i,j) += this->operator()(i,k)*other.operator()(k,j);
            return mult;
        }
        matrix operator  /(T d)
        {
            assert(d != (T)(0));
            matrix m(this->N, this->M);
            for(std::size_t i = 0; i < m.N; ++i) for(std::size_t j = 0; j < m.M; ++j)
                m.operator()(i,j) = this->operator()(i,j)/d;
            return m;
        }
        void   operator +=(matrix other)
        {
            assert(this->N == other.N && this->M == other.M);
            for (std::size_t i = 0; i < this->N; ++i) for (std::size_t j = 0; j < this->M; ++j)
                this->operator()(i,j)+=other.operator()(i,j);
            return;
        }
        void   operator -=(matrix other)
        {
            assert(this->N == other.N && this->M == other.M);
            for (std::size_t i = 0; i < this->N; ++i) for (std::size_t j = 0; j < this->M; ++j)
                this->operator()(i,j)-=other.operator()(i,j);
            return;
        }
        void   operator *=(matrix other)
        {
            assert(this->M == other.N || (other.size() == 1));
            if (other.size() == 1)
            {
                matrix mult(this->N, this->M);
                for(std::size_t i = 0; i < mult.N; ++i)
                    for(std::size_t j = 0; j < mult.M; ++j)
                        mult.operator()(i,j) = other.operator()(0,0)*this->operator()(i,j);
                *this = mult;
                return;
            }
            matrix mult(this->N, other.M);
            for(std::size_t i = 0; i < mult.N; ++i)
                for(std::size_t j = 0; j < mult.M; ++j)
                    for(std::size_t k = 0; k < this->M; ++k)
                        mult.operator()(i,j) += this->operator()(i,k)*other.operator()(k,j);
            *this = mult;
            return;
        }
        void   operator /=(T d)
        {
            assert(d != (T)(0));
            matrix<T>m(this->M,this->M);
            for (std::size_t i = 0; i < m.N; ++i) m(i,i) = 1/d;
            this->operator*=(m);
        }
        matrix transpose()
        {
            matrix m(this->M,this->N);
            for(std::size_t i = 0; i < this->M; ++i) for(std::size_t j = 0; j < this->N; ++j)
                m(i,j) = this->operator()(j,i);
            return m;
        }
        matrix column(std::size_t j)
        {
            assert(j >= 0 && j < this->M);
            matrix sub(this->N,1);
            for(std::size_t i = 0; i < this->N; ++i)
                sub.operator()(i,0) = this->operator()(i,j);
            return sub;
        }
        matrix row   (std::size_t i)
        {
            assert(i >= 0 && i < this->N);
            matrix sub(1,this->M);
            for(std::size_t j = 0; j < this->M; ++j) {std::cout << i << ' ' << j << std::endl; sub.operator()(0,j) = this->operator()(i,j);}
            return sub;
        }
        T conv()
        {
            assert(this->N == 1 && this->M == 1);
            return this->operator()(0,0);
        }
        void print()
        {
            for (std::size_t i = 0; i < this->N; ++i)
            {
                for (std::size_t j = 0; j < this->M; ++j) std::cout << this->operator()(i,j) << ' ';
                std::cout << std::endl;
            }
            return;
        }
    };

    template<typename T>
    matrix<T> diag(T d, std::size_t N)
    {
        matrix<T>m(N,N);
        for (std::size_t i = 0; i < m.size(); ++i) m(i,i) = d;
        return m;
    }
    template<typename T>
    matrix<T> Id(std::size_t N)
    {
        return diag<T>((T)(1), N);
    }
    template<typename T>
    matrix<T> E(std::size_t i, std::size_t j, std::size_t N)
    {
        assert(i >= 0 && j >= 0 && i <= N && j <= N);
        matrix<T> m(N); m(i, j) = 1;
        return m;
    }

    template<typename T>
    class spatial_vector: public matrix<T>
    {
        public:
        T x, y, z;

        private:
        void correct_xyz_parameters()
        {
            this->x = this->operator()(0);
            this->y = this->operator()(1);
            this->z = this->operator()(2);
            return;
        }

        public:

        spatial_vector(): matrix<T>(3,1)
        {
            this->operator()(0) = 0; this->x = 0;
            this->operator()(1) = 0; this->y = 0;
            this->operator()(2) = 0; this->z = 0;
        }
        spatial_vector(T x_, T y_, T z_): matrix<T>(3,1)
        {
            this->operator()(0) = x_; this->x = x_;
            this->operator()(1) = y_; this->y = y_;
            this->operator()(2) = z_; this->z = z_;
        }
        spatial_vector(matrix<T> m): matrix<T>(3,1)
        {
            assert(m.size() == 3);
            this->operator()(0) = m.operator()(0); this->x = m.operator()(0);
            this->operator()(1) = m.operator()(1); this->y = m.operator()(1);
            this->operator()(2) = m.operator()(2); this->x = m.operator()(2);
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
            this->correct_xyz_parameters();
            return;
        }
        void operator *=(T d)
        {
            this->operator()(0) *= d;
            this->operator()(1) *= d;
            this->operator()(2) *= d;
            this->correct_xyz_parameters();
        }
        void operator /=(T d)
        {
            this->operator()(0) /= d;
            this->operator()(1) /= d;
            this->operator()(2) /= d;
            this->correct_xyz_parameters();
        }
        T    operator  *(spatial_vector other)
        {
            return ((matrix(this->data,3,1).transpose()).operator*(other)).conv();
        }
        T norm()
        {
            return sqrt(this->operator*(*this));
        }
        spatial_vector operator^(spatial_vector other)
        {
            return spatial_vector(
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
            this->correct_xyz_parameters();
        }
        void print()
        {
            std::cout << '(' << this->operator()(0) << ", " << 
                                this->operator()(1) << ", " <<
                                this->operator()(2) << ')' << std::endl;
        }
    };

    class quaternion
    {
        public:

        double scalar;
        spatial_vector<double> v;
    };
};