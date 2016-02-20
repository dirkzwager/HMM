#include <iostream>

const float PI = 3.14;

using namespace std;

template <class T>
class Matrix {
public:
    Matrix() { mtrx=NULL; rows=0; cols=0; }
    ~Matrix();
    Matrix(int ncols, int nrows);
    Matrix(float anglex, float angley, float anglez);
    Matrix<T> multMatrix(Matrix m2, int nrowsB);
    Matrix<T> multMatrix(Matrix m2);
    Matrix<T> squareMatrix(int power);                
    void fillMatrix(T* vals, int N);
    void printMatrix();
    T* flatMatrix();
    Matrix<T> rotationXMatrix3D(float angle);    
    Matrix<T> rotationYMatrix3D(float angle);
    Matrix<T> rotationZMatrix3D(float angle);
    Matrix<T> rotationMatrix3D(float anglex, float angley, float anglez);

    T cellVal(const int i, const int j) { return mtrx[i][j]; }
    void setVal(const int i, const int j, T val) { mtrx[i][j] = val; }
    Matrix<T>  col(const int c);
    Matrix<T>  row(const int r);
    T  		   sumCol(const int c);
    T          sumRow(const int r);
    
    Matrix<T>& operator =    (const Matrix<T>& m2);
    Matrix<T>  operator *    (const Matrix<T>& m2);
    Matrix<T>& operator *=   (Matrix<T>& m2);
    Matrix<T>  operator ^    (const int pow);
    Matrix<T>& operator ^=   (const int pow);
private:        
    T** mtrx;
    int rows;
    int cols;
};

template <class T>
Matrix<T>::Matrix(int ncols, int nrows) {
    rows = nrows;
    cols = ncols;
    mtrx = new T*[ncols];
    
    for (int i=0; i<cols; i++) 
        mtrx[i] = new T[rows];   
        
    for (int i=0; i<cols; i++)         
        for (int j=0; j<rows; j++) 
            mtrx[i][j] = (T)0.0;       
}

template <class T>
Matrix<T>::~Matrix() {
	/*for (int i=0; i<cols; i++)
		delete[] mtrx[i];
	delete[] mtrx;*/
}

template <class T>
void Matrix<T>::fillMatrix(T* vals, int N) {
    if (rows*cols == N) {
        int r = 0;
        int c = 0;
        
        for (int i=0; i<N; i++) {
            mtrx[c][r] = vals[i];
            if (r<rows-1)
                r++;
            else {
                r=0;
                c++;
            }
        }
    }
}

//mult mtrx with m2, works only for rows == cols.
template <class T>
Matrix<T> Matrix<T>::multMatrix(Matrix m2, int nrowsB) {
    Matrix<T> matrix(cols, nrowsB);// = new int*[cols];
            
    for (int i=0; i<nrowsB; i++) 
        for (int j=0; j<cols; j++) 
            for (int k=0; k<cols; k++) 
                matrix.mtrx[i][j] += mtrx[i][k] * m2.mtrx[k][j];            
    return matrix;
}

//mult mtrx with m2
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m2) {
	if (mtrx!=NULL) {
		for (int i=0; i<cols; i++)
			delete[] mtrx[i];
		delete[] mtrx;
	}
	mtrx = m2.mtrx;
	cols = m2.cols;
	rows = m2.rows;
	return *this;
}

//mult mtrx with m2
template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m2) {
	if (cols == m2.rows && rows == m2.cols) {
	    Matrix<T> m(m2.cols, rows);
    	for (int r1=0; r1<rows; r1++) 
    		for (int c2=0; c2<m2.cols; c2++) 
    			for (int c1=0; c1<cols; c1++) 
    				m.mtrx[c2][r1] += mtrx[c1][r1] * m2.mtrx[c2][c1];
    	return m;
    }
    else printf("\n\nError, incorrect matrix sizes..\n\n");
}

//mult mtrx with m2
template <class T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T>& m2) {
	if (cols == m2.rows && rows == m2.cols) {
	    Matrix<T> m(m2.cols, rows);
    	for (int r1=0; r1<rows; r1++) 
    		for (int c2=0; c2<m2.cols; c2++) 
    			for (int c1=0; c1<cols; c1++) 
    				m.mtrx[c2][r1] += mtrx[c1][r1] * m2.mtrx[c2][c1];
    	*this = m;
    	return *this;
    }
    else printf("\n\nError, incorrect matrix sizes..\n\n");
}

//mult mtrx with m2
template <class T>
Matrix<T>& Matrix<T>::operator^=(const int power) {
		for (int i=0; i<power-1; i++)
			*this*=(*this);
		return *this;
}

//mult mtrx with m2
template <class T>
Matrix<T> Matrix<T>::multMatrix(Matrix m2) {
	if (cols == m2.rows && rows == m2.cols) {
	    Matrix<T> m(m2.cols, rows);
    	for (int r1=0; r1<rows; r1++) 
    		for (int c2=0; c2<m2.cols; c2++) 
    			for (int c1=0; c1<cols; c1++) 
    				m.mtrx[c2][r1] += mtrx[c1][r1] * m2.mtrx[c2][c1];
    	return m;
    }
    else printf("\n\nError, incorrect matrix sizes..\n\n");
}

template <class T>
T* Matrix<T>::flatMatrix() {
    T* flam = new T[rows*cols];
    int idx = 0;
    for (int c=0; c<cols; c++)
        for (int r=0; r<rows; r++) {   
            flam[idx] = mtrx[c][r] ;       
            idx++;
        }    
    return flam;
}

template <class T>
Matrix<T> Matrix<T>::col(const int c) {
	Matrix<T> m(1, rows);
	for (int r=0; r<rows; r++)
		m.setVal(0,r,mtrx[c][r]);
	return m;
}

template <class T>
T Matrix<T>::sumCol(const int c) {
	T sum = (T)0;
	for (int r=0; r<rows; r++)
		sum += mtrx[c][r];
	return sum;
}

template <class T>
Matrix<T> Matrix<T>::row(const int r) {
	Matrix<T> m(cols, 1);
	for (int c=0; c<cols; c++)
		m.setVal(c,0,mtrx[c][r]);
	return m;
}

template <class T>
T Matrix<T>::sumRow(const int r) {
	T sum = (T)0;
	for (int c=0; c<cols; c++)
		sum += mtrx[c][r];
	return sum;
}

template <class T>
Matrix<T> Matrix<T>::rotationXMatrix3D(float angle) {
    Matrix<T> rxmatrix(4,4);
    float rad = (PI/180.0)*angle;
    float* vals = new float[16];
        
    vals[0] = 1.0;  vals[4] = 0.0;         vals[8] = 0.0;             vals[12] = 0.0;
    vals[1] = 0.0;  vals[5] = cos(rad);    vals[9] = -1.0*sin(rad);   vals[13] = 0.0;
    vals[2] = 0.0;  vals[6] = sin(rad);    vals[10] = cos(rad);       vals[14] = 0.0;
    vals[3] = 0.0;  vals[7] = 0.0;         vals[11] = 0.0;            vals[15] = 1.0;        
    rxmatrix.fillMatrix(vals, 16);
    return rxmatrix;
}

template <class T> 
Matrix<T> Matrix<T>::rotationYMatrix3D(float angle) {
    Matrix<T> rymatrix(4,4);
    float rad = (PI/180)*angle;
    float* vals = new float[16];
        
    vals[0] = cos(rad);      vals[4] = 0.0;    vals[8] = sin(rad);     vals[12] = 0.0;
    vals[1] = 0.0;           vals[5] = 1.0;    vals[9] = 0.0;          vals[13] = 0.0;
    vals[2] = -1.0*sin(rad); vals[6] = 0.0;    vals[10] = cos(rad);    vals[14] = 0.0;
    vals[3] = 0.0;           vals[7] = 0.0;    vals[11] = 0.0;         vals[15] = 1.0;        
    rymatrix.fillMatrix(vals, 16);
    return rymatrix;
}

template <class T>
Matrix<T> Matrix<T>::rotationZMatrix3D(float angle) {
    Matrix<T> rzmatrix(4,4);
    float rad = (PI/180)*angle;
    float* vals = new float[16];
        
    vals[0] = cos(rad); vals[4] = -1.0*sin(rad); vals[8] = 0.0;    vals[12] = 0.0;
    vals[1] = sin(rad); vals[5] = cos(rad);      vals[9] = 0.0;    vals[13] = 0.0;
    vals[2] = 0.0;      vals[6] = 0.0;           vals[10] =1.0;    vals[14] = 0.0;
    vals[3] = 0.0;      vals[7] = 0.0;           vals[11] = 0.0;   vals[15] = 1.0;        
    rzmatrix.fillMatrix(vals, 16);
    return rzmatrix;
}

template <class T>
Matrix<T>::Matrix(float anglex, float angley, float anglez) {
    Matrix rmtrx = rotationXMatrix3D(anglex).multMatrix(rotationYMatrix3D(angley),4);
    rows = 4;
    cols = 4;
    mtrx = new T*[cols];
    for (int i=0; i<cols; i++) 
        mtrx[i] = new T[rows];
    for (int i=0; i<cols; i++)         
        for (int j=0; j<rows; j++) 
            mtrx[i][j] = rmtrx.mtrx[i][j]; 
}

template <class T>
void Matrix<T>::printMatrix() {
	for (int i=0; i<rows; i++) {
		cout << "| ";
		for (int j=0; j<cols; j++) {
			if (j==cols-1)
				cout << mtrx[j][i] << " |";
			else
				cout << mtrx[j][i] << "\t"; 
		}
		cout << endl;
	}
}
