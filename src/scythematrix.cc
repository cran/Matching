/* Edited by Jasjeet S. Sekhon <jasjeet_sekhon@harvard.edu> */
/* HTTP://jsekhon.fas.harvard.edu/                          */
/* Feb 14,  2005                                            */
//
// Memeber function definitions for the Scythe_Double_Matrix.h
// header file.  These functions make up the Matrix class as used
// in the Scythe project.  
//
// Scythe C++ Library
// Copyright (C) 2000 Kevin M. Quinn and Andrew D. Martin
//
// This code written by:
//
// Kevin Quinn
// Assistant Professor
// Dept. of Political Science and 
// Center for Statistics and the Social Sciences
// Box 354322
// University of Washington
// Seattle, WA  98195-4322
// quinn@stat.washington.edu
//
// Andrew D. Martin
// Assistant Professor
// Dept. of Political Science
// Campus Box 1063
// Washington University
// St. Louis, MO 63130
// admartin@artsci.wustl.edu
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef SCYTHE_DOUBLE_MATRIX_CC
#define SCYTHE_DOUBLE_MATRIX_CC

#include "scythematrix.h"

/*! \class Matrix
 *  \brief The Matrix Class that is the cornerstone of the Scythe 
 *	Statistical Library
 
 *  The Matrix class contains functions that modify the Matrix object 
 *	in a variety of ways.
 */
// Avoid NameSpace Pollution
namespace SCYTHE {
    using namespace std;
    
/********************     BASIC FUNCTIONS     *************************/
/**********************************************************************/
    
/*!
 * \brief CONSTRUCTOR:  Matrix  - Creates Matrix obj. with specificed 
 * rows and columns, sets each element = 0
 * \param rows An integer that reflects the number of rows.
 * \param cols An integer that reflects the numbef of columns.
 */
    Matrix::Matrix (const int& rows, const int& cols)
    {
	if (rows < 1 || cols < 1) {
	    cerr << "ERROR 0001: Improper row or column dimension in "
		 << "Matrix constructor"
		 << endl;
	    exit (1);
	}
	rowsize = rows;    // assign Matrix rowsize
	colsize = cols;    // assign Matrix colsize 
	size = rows * cols;    // assign Matrix size
	data = new double[size];
	
	for (int i = 0; i < size; ++i) {
	    data[i] = 0.0;
	}
    }
    
/*!
 * \brief CONSTRUCTOR:  Matrix - Creates Matrix object filled with data 
 * from array
 * \param inputarray An array containing the data for the Matrix in 
 * row major order.
 * \param rows An integer reflecting the number of rows.
 * \param cols An integer that reflects the numbef of columns.
 */
    Matrix::Matrix (const double *inputarray, const int& rows, const int& cols)
    {
	if (rows < 1 || cols < 1) {
	    cerr << "ERROR 0002: Improper row or column dimension "
		 << "in Matrix constructor"
		 << endl;
	    exit (2);
	}
	rowsize = rows;    // assign Matrix rowsize
	colsize = cols;    // assign Matrix colsize
	size = rows * cols;    // assign Matrix size
	data = new double[size];
	
	for (int i = 0; i < size; ++i) {
	    data[i] = inputarray[i];
	}
    }
    
/*! 
 * \brief CONSTRUCTOR:  Matrix - Creates Matrix object from old Matrix
 * \param old_Matrix A Matrix object that contains the data from 
 * another Matrix.
 */
    Matrix::Matrix (const Matrix & old_Matrix)
    {
	rowsize = old_Matrix.rowsize;  // assign Matrix rowsize
	colsize = old_Matrix.colsize;  // assign Matrix colsize
	size = old_Matrix.size;  // assign Matrix size
	data = new double[size];
	
	for (int i = 0; i < size; ++i) {
	    data[i] = old_Matrix.data[i];
	}
    }
    
/*!
 * \brief OPERATOR:  Matrix operator '=' - Allows for the copying of a 
 * Matrix
 * \param B A Matrix from which the data will be copied.
 * \return A Matrix object.
 */
    
    Matrix & 
    Matrix::operator = (const Matrix & B)
    {
	rowsize = B.rowsize;
	colsize = B.colsize;
	size = B.size;
	
	delete[]data;
	data = new double[size];
	for (int i = 0; i < size; ++i)
	    data[i] = B.data[i];
	
	return *this;
    }
    
    
    
    
/*!
 * \brief OPERATOR: Matrix operator () - Retrieves all Matrix elements 
 * in row \a i.  Please Note: Indexing starts at 0.
 * \param i a constant integer referring to the number of the 
 * row to be retrived.
 * \param a a constant all_elements 
 *  \return Matrix all elements in row \a i.
 */
    
    Matrix 
    Matrix::operator () (const int& i, const all_elements& a)
    {
	if (i >= rowsize || i < 0) {
	    cerr << "ERROR 0005: Index out of range in () operator" << endl;
	    exit (5);
	}
	
	int newrowsize = 1;
	int newcolsize = colsize;
	double *newdata = new double[newcolsize];
	
	for (int j = 0; j < newcolsize; ++j) {
	    newdata[j] = data[j + i*colsize];
	}
	
	Matrix temp = Matrix (newdata, newrowsize, newcolsize);
	delete[]newdata;
	return temp;
	
    }
    
    
    
/*! 
 * \brief OPERATOR: Matrix operator () - Retrieves all Matrix elements 
 * in column \a j. Please Note: Indexing starts at 0.
 * \param _ a constant of type all_elements.
 * \param j a constant integer referring to the number of the 
 * column to be retrived.
 *  \return Matrix all elements in column \a j.
 */
    
    Matrix 
    Matrix::operator () (const all_elements& a, const int& j)
    {
	if (j >= colsize || j < 0) {
	    cerr << "ERROR 0008: Index out of range in () operator" << endl;
	    exit (8);
	}
	
	int newrowsize = rowsize;
	int newcolsize = 1;
	double *newdata = new double[newrowsize];
	
	for (int i = 0; i < newrowsize; ++i) {
	    newdata[i] = data[j + i*colsize];
	}
	
	Matrix temp = Matrix (newdata, newrowsize, newcolsize);
	delete[]newdata;
	return temp;
    }
    
    
    
/*!
 * \brief OPERATOR: Matrix operator () - Retrieves Matrix elements 
 * in row \a i.  Please Note: Indexing starts at 0.
 * \param i a constant integer referring to the number of the 
 * row to be retrived.
 * \param J a constant reference to the Matrix from which the 
 * data will be extracted.
 *  \return Matrix elements in row \a i.
 */
    
    Matrix
    Matrix::operator () (const int& i, const Matrix& J)
    {
	
	if (i >= rowsize || i < 0) {
	    cerr << "ERROR 0005: Index out of range in () operator" << endl;
	    exit (5);
	}
	
	if (J.colsize != 1 && J.rowsize != 1) {
	    cerr << "ERROR 0006: Either rows or cols of J != 1 in () operator" 
		 << endl;
	    exit (6);
	}
	
	int newrowsize = 1;
	int newcolsize = J.size;
	double *newdata = new double[newcolsize];
	
	for (int j = 0; j < newcolsize; ++j) {
	    int index = static_cast < int >(J.data[j]);
	    if (index >= colsize || index < 0) {
		cerr << "ERROR 0007: Index out of range in () operator" << endl;
		exit (7);
	    }
	    index = index + i * colsize;
	    newdata[j] = data[index];
	}
	
	Matrix temp = Matrix (newdata, newrowsize, newcolsize);
	delete[]newdata;
	return temp;
    }
    
/*! 
 * \brief OPERATOR: Matrix operator () - Retrieves Matrix elements 
 * in column \a j. Please Note: Indexing starts at 0.
 * \param I a constant reference to the Matrix from which the data 
 * will be extracted.
 * \param j a constant integer referring to the number of the 
 * column to be retrived.
 *  \return Matrix elements in column \a j.
 */
    Matrix
    Matrix::operator () (const Matrix& I, const int& j)
    {
	
	if (j >= colsize || j < 0) {
	    cerr << "ERROR 0008: Index out of range in () operator" << endl;
	    exit (8);
	}
	
	if (I.colsize != 1 && I.rowsize != 1) {
	    cerr << "ERROR 0009: Either rows or cols of I != 1 in () operator" << endl;
	    exit (9);
	}
	
	int newrowsize = I.size;
	int newcolsize = 1;
	double *newdata = new double[newrowsize];
	
	for (int i = 0; i < newrowsize; ++i) {
	    int index = static_cast < int >(I.data[i]);
	    if (index >= rowsize || index < 0) {
		cerr << "ERROR 0010: Index out of range in () operator" << endl;
		exit (10);
	    }
	    index = j + index * colsize;
	    newdata[i] = data[index];
	}
	
	Matrix temp = Matrix (newdata, newrowsize, newcolsize);
	delete[]newdata;
	return temp;
    }
    
    
    
/*! 
 * \brief OPERATOR: Matrix operator () - Extracts submatrix 
 * from existing matrix
 
 * OPERATOR: Matrix operator () - Extracts submatrix from existing matrix.  
 * Get elements \a i,\a j from a Matrix where \a i in \a I and \a j 
 * in \a J.  Indexing starts at 0.
 * \param I a constant reference to the Matrix from which the data 
 * will be extracted.
 * \param J a constant reference to another Matrix from which the 
 * data will be extracted.
 * \return a new Matrix created from the selected data from the 
 * previous two.
 */
    Matrix 
    Matrix::operator () (const Matrix& I, const Matrix& J){
	if (I.colsize != 1 && I.rowsize != 1) {
	    cerr << "ERROR 0011: Either Rows or Cols of I != 1 in () operator" 
		 << endl;
	    exit (11);
	}
	if (J.colsize != 1 && J.rowsize != 1) {
	    cerr << "ERROR 0012: either rows or cols of J != 1 in () operator"
		 << endl;
	    exit (12);
	}
	if (I.size > rowsize){
	    cerr << "ERROR 0013: size(I) > rowsize of Matrix in Matrix operator ()"
		 << endl;
	    exit(13);
	}
	if (J.size > colsize){
	    cerr << "ERROR 0014: size(J) > colsize of Matrix in Matrix operator ()"
		 << endl;
	    exit(14);
	}
	
	int place = 0;
	int indexi, indexj;
	double *newdata = new double[I.size * J.size];
	for (int i = 0; i < I.size; i++) {
	    for (int j = 0; j < J.size; j++) {
		indexi = static_cast < int > (I.data[i]);
		indexj = static_cast < int > (J.data[j]);
		if (indexi >= rowsize || indexi < 0) {
		    cerr << "ERROR 0016: Row index out of range in () operator" 
			 << endl;
		    exit (16);
		}
		if (indexj >= colsize || indexj < 0) {
		    cerr << "ERROR 0017: Column index out of range in () operator"
			 << endl;
		    exit (17);
		}
		newdata[place] = data[indexi * colsize + indexj];
		place++;
	    }
	}
	
	Matrix temp = Matrix (newdata, I.size, J.size);
	delete[]newdata;
	return temp;
    }
    
    
    
/*!
 * \brief Prints the Matrix to the screen
 * \param width a constant integer reflecting the screen width 
 * (in characters) for each element of the Matrix.  This is used 
 * in the C++ \e setw() function.
 * \param prec a constant integer reflecting the decimal precision 
 * of each element in the Matrix.  This is used in the C++ \e 
 * setprecision() function.
 * \return void
 */
    void
    Matrix::print (const int width, const int prec)
    {
	
	int count = 0;
	
	for (int i = 0; i < rowsize; ++i) {
	    if (i > 0) {
		cout << endl;
	    }
	    for (int j = 0; j < colsize; ++j) {
		cout << setw (width) << setprecision (prec) 
		     << data[i * colsize + j] << " ";
		++count;
	    }
	}
	cout << endl << endl;
    }
    
/********************   MORE ADVANCED FUNCTIONS  **********************/
/**********************************************************************/
    
//  FUNCTION: c  - concatenates a sequence of doubles into a Matrix
/*!
 * \brief Concatenates a sequence of doubles into a Matrix.
 * 
 * Concatenates a sequence of doubles into a Matrix.
 * \param a first double to be concatenated.
 * \param b second double to be concatenated.
 * \param ... other doubles (number determined by user (up to 26)).
 * \return A Matrix object, the column vector formed by concatenating the 
 * input doubles.
 */
    Matrix c (const double& a, const double& b){
	double *newdata = new double[2];
	newdata[0] = a;
	newdata[1] = b;
	
	Matrix temp = Matrix (newdata, 2, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c){
	double *newdata = new double[3];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	
	Matrix temp = Matrix (newdata, 3, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d){ 
	double *newdata = new double[4];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	
	Matrix temp = Matrix (newdata, 4, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e){ 
	double *newdata = new double[5];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	
	Matrix temp = Matrix (newdata, 5, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f){ 
	double *newdata = new double[6];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	
	Matrix temp = Matrix (newdata, 6, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g){ 
	double *newdata = new double[7];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	
	Matrix temp = Matrix (newdata, 7, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h){ 
	double *newdata = new double[8];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	
	Matrix temp = Matrix (newdata, 8, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i){ 
	double *newdata = new double[9];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	
	Matrix temp = Matrix (newdata, 9, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j){ 
	double *newdata = new double[10];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	
	Matrix temp = Matrix (newdata, 10, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k){ 
	double *newdata = new double[11];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	
	Matrix temp = Matrix (newdata, 11, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l){ 
	double *newdata = new double[12];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	
	Matrix temp = Matrix (newdata, 12, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m){ 
	double *newdata = new double[13];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	
	Matrix temp = Matrix (newdata, 13, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n){ 
	double *newdata = new double[14];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	
	Matrix temp = Matrix (newdata, 14, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o){ 
	double *newdata = new double[15];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	
	Matrix temp = Matrix (newdata, 15, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p){ 
	double *newdata = new double[16];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	
	Matrix temp = Matrix (newdata, 16, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q){ 
	double *newdata = new double[17];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	
	Matrix temp = Matrix (newdata, 17, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r){ 
	double *newdata = new double[18];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	
	Matrix temp = Matrix (newdata, 18, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s){ 
	double *newdata = new double[19];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	
	Matrix temp = Matrix (newdata, 19, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t){ 
	double *newdata = new double[20];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	
	Matrix temp = Matrix (newdata, 20, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u){ 
	double *newdata = new double[21];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	
	Matrix temp = Matrix (newdata, 21, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v){ 
	double *newdata = new double[22];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	newdata[21] = v;
	
	Matrix temp = Matrix (newdata, 22, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w){ 
	double *newdata = new double[23];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	newdata[21] = v;
	newdata[22] = w;
	
	Matrix temp = Matrix (newdata, 23, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x){ 
	double *newdata = new double[24];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	newdata[21] = v;
	newdata[22] = w;
	newdata[23] = x;
	
	Matrix temp = Matrix (newdata, 24, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x,
	      const double& y){ 
	double *newdata = new double[25];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	newdata[21] = v;
	newdata[22] = w;
	newdata[23] = x;
	newdata[24] = y;
	
	Matrix temp = Matrix (newdata, 25, 1);
	delete[]newdata;
	return temp;
    }
    
    Matrix c (const double& a, const double& b, const double& c, 
	      const double& d, const double& e, const double& f,
	      const double& g, const double& h, const double& i,
	      const double& j, const double& k, const double& l,
	      const double& m, const double& n, const double& o,
	      const double& p, const double& q, const double& r,
	      const double& s, const double& t, const double& u,
	      const double& v, const double& w, const double& x,
	      const double& y, const double& z){ 
	double *newdata = new double[26];
	newdata[0] = a;
	newdata[1] = b;
	newdata[2] = c;
	newdata[3] = d;
	newdata[4] = e;
	newdata[5] = f;
	newdata[6] = g;
	newdata[7] = h;
	newdata[8] = i;
	newdata[9] = j;
	newdata[10] = k;
	newdata[11] = l;
	newdata[12] = m;
	newdata[13] = n;
	newdata[14] = o;
	newdata[15] = p;
	newdata[16] = q;
	newdata[17] = r;
	newdata[18] = s;
	newdata[19] = t;
	newdata[20] = u;
	newdata[21] = v;
	newdata[22] = w;
	newdata[23] = x;
	newdata[24] = y;
	newdata[25] = z;
	
	Matrix temp = Matrix (newdata, 26, 1);
	delete[]newdata;
	return temp;
    }
    
    
    
//  FUNCTION: Transpose  - computes the transpose of a Matrix
/*!
 * \brief Computes the transpose of the given Matrix.
 * 
 * Computes the transpose of the given Matrix.
 * \param old_matrix a constant reference to a Matrix object.  
 * This Matrix will be transposed.
 * \return A Matrix object, the transpose of the input Matrix object.
 */
    Matrix t (const Matrix & old_matrix)
    {
	int newrowsize = old_matrix.colsize;
	int newcolsize = old_matrix.rowsize;
	double *newdata = new double[old_matrix.size];
	for (int i = 0; i < newcolsize; ++i) {
	    for (int j = 0; j < newrowsize; ++j) {
		newdata[i + newcolsize * j] = old_matrix.data[j + newrowsize * i];
	    }
	}
	Matrix temp = Matrix (newdata, newrowsize, newcolsize);
	delete[]newdata;
	return temp;
    }
    
/*!
 * \brief Creates a Matrix of Ones
 *
 * Creates a Matrix filled with Ones, given a specified size.
 * \param rows a constant int reflecting the number of rows in the Matrix.
 * \param cols a constant int reflecting the number of columns in the Matrix.
 * \return a new Matrix filled with 1's.
 */
    Matrix ones (const int& rows, const int& cols)
    {
	if (rows < 1 || cols < 1) {
	    cerr << "Error 0018: improper row or column dimension in ones()"
		 << endl;
	    exit (18);
	}
	double *newdata = new double[rows * cols];
	int size = rows * cols;
	for (int i = 0; i < size; ++i) {
	    newdata[i] = 1.0;
	}
	Matrix temp = Matrix (newdata, rows, cols);
	delete[]newdata;
	return temp;
    }
    
/*!
 * \brief Creates a Matrix of Zeros
 *
 * Creates a Matrix filled with Zeros, given a specified size.
 * \param rows a constant int reflecting the number of rows in the Matrix.
 * \param cols a constant int reflecting the number of columns in the Matrix.
 * \return a new Matrix filled with 1's.
 */
    Matrix zeros (const int& rows, const int& cols)
    {
	if (rows < 1 || cols < 1) {
	    cerr << "Error 0018: improper row or column dimension in ones()"
		 << endl;
	    exit (18);
	}
	double *newdata = new double[rows * cols];
	int size = rows * cols;
	for (int i = 0; i < size; ++i) {
	    newdata[i] = 0.0;
	}
	Matrix temp = Matrix (newdata, rows, cols);
	delete[]newdata;
	return temp;
    } // end of zeros
    
    
    
//  FUNCTION: Eye - creates an Identity Matrix of size k x k
/*! 
 * \brief Creates an Identity Matrix
 *
 * Creates an Identity Matrix of size \a k \a x \a k.
 * \param k a constant integer reflecting the length and width 
 * of the identity matrix.
 * \return the Identity Matrix
 */
    Matrix eye (const int& k)
    {
	double *newdata = new double[k * k];
	double hold;
	for (int i = 0; i < k; ++i) {
	    for (int j = 0; j < k; ++j) {
		if (i == j)
		    hold = 1.0;
		else
		    hold = 0.0;
		newdata[k * i + j] = hold;
	    }
	}
	Matrix temp = Matrix (newdata, k, k);
	delete[]newdata;
	return temp;
    }
    
/*! 
 * \brief Creates a Vector-additive sequence Matrix
 *
 * Creates a Vector-additive sequence Matrix of (\a size x 1)
 * \param start a constant double reflecting the start value of 
 * the first element in the vector.
 * \param incr a double constant reflecting the incremental step 
 * value between each matrix element.
 * \param size a constant integer reflecting the size of the vector.
 * \return a new Matrix (vector).
 */
    Matrix seqa (const double& start, const double& incr, const int& size)
    {
	double *newdata = new double[size];
	double val = start;
	for (int i = 0; i < size; ++i) {
	    newdata[i] = val;
	    val += incr;
	}
	Matrix temp = Matrix (newdata, size, 1);
	delete[]newdata;
	return temp;
    }
    
    
/*! 
 * \brief Sorts all elements of a Matrix (not column by column) using 
 * shellsort
 *
 * Sorts all elements of a Matrix (not column by column) using shellsort
 * \param A the Matrix to be sorted.
 * \return a new Matrix the same size as the original in 
 * which all elements have been sorted.
 */
    Matrix sort(const Matrix& A){
	int i, j, h;
	double v;
	
	double *newdata = new double[A.size];
	for (i = 0; i<A.size; ++i)
	    newdata[i] = A.data[i];
	for (h = 1; h <= A.size/9; h = 3*h+1);
	for (; h > 0; h /= 3)
	    for (i = h+1; i <= A.size; i += 1){
		v = newdata[i-1]; 
		j = i;
		while (j>h && newdata[j-h-1] > v)
		{newdata[j-1] = newdata[j-h-1]; j -= h;}
		newdata[j-1] = v;
	    }
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;     
    }
    
    
    
/*! 
 * \brief Sorts all columns of a Matrix using shellsort
 *
 * Sorts all columns of a Matrix using shellsort
 * \param A the Matrix to be sorted.
 * \return a new Matrix the same size as the original in 
 * which all elements have been sorted.
 */
    Matrix sortc(const Matrix& A){
	int i, j, h;
	double v;
	
	double *newdata = new double[A.size];
	for (i = 0; i<A.size; ++i)
	    newdata[i] = A.data[i];
	
	for (int colindex=0; colindex<A.colsize; ++colindex){
	    
	    for (h = 1; h <= A.rowsize/9; h = 3*h+1);
	    for (; h > 0; h /= 3)
		for (i = h+1; i <= A.rowsize; i += 1){
		    v = newdata[(i-1)*A.colsize + colindex]; 
		    j = i;
		    while (j>h && newdata[(j-h-1)*A.colsize + colindex] > v){
			newdata[(j-1)*A.colsize + colindex] = 
			    newdata[(j-h-1)*A.colsize + colindex]; 
			j -= h;
		    }
		    newdata[(j-1)*A.colsize + colindex] = v;
		}
	    
	}
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;     
    }
    
    
/*! 
 * \brief Cholesky Decomposition
 *
 * Cholesky Decomposition of a Symmetric Positive Definite Matrix. 
 * Given an input Matrix \a A computes \a L where \a L*L' \a = \a A.  
 * This function returns the lower triangular matrix \a L.  
 * NOTE:  This function does not check for symmetry.
 * \param a A constant reference to a Matrix.
 * \return the new Cholesky Decomposition of the input Matrix 
 * (the lower triangular Matrix \a L).
 */
    Matrix cholesky (const Matrix & A){
	
	if (A.rowsize == A.colsize){
	    register double *newdata = new double[A.rowsize * A.colsize];
	    for (int i = 0; i < A.rowsize; ++i){
		for (int j = i; j < A.colsize; ++j){
		    register double h = A.data[A.colsize * i + j];
		    for (int k = 0; k <= i - 1; ++k){
			h -= newdata[A.colsize * i + k] * newdata[A.colsize * j + k];
		    }
		    if (i == j){
			if (h <= 0){
			    cerr << "ERROR 0019: input Matrix not positive definite in SCYTHE::cholesky()\n"
				 << "Cholesky decomposition failed. Exiting program. "
				 << endl;
			    exit (19);
			}
			newdata[A.colsize * i + i] = ::sqrt (h);
		    } else {
			newdata[A.colsize * j + i] = (1.0 / newdata[A.colsize * i + i])*h;
			newdata[A.colsize * i + j] = 0.0;
		    }
		}
	    }
	    
	    Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	    delete[]newdata;
	    return temp;
	} else {
	    cerr << "ERROR 0020: input Matrix not square in SCYTHE::cholesky() \n" 
		 << "Cholesky decomposition failed. Exiting program. "
		 << endl;
	    exit (20);
	}
    }
    
/*! 
 * \brief Solves \a A \a x = \a b for x via back-substitution 
 * using Cholesky Decomposition
 *
 * Solves \a A \a x = \a b for x via back-substitution using 
 * Cholesky Decomposition.  \a A must be symmetric and positive definite.
 * \remarks The solution is trivial because of the following identity.
 * A*x = (M*M')*x = M*(M'*x) = b.
 * Function works by successively solving two triangular systems:
 * M*y = b for y
 * and then
 * M'*x = y for x.
 * \param A a constant reference to the Matrix \e A in the equality 
 * \a A \a x = \a b.
 * \param b a constant reference to the Matrix \e b in the equality
 * \a A \a x = \a b.
 * \return a new Matrix (vector) \a x from the equality
 * \a A \a x = \a b.
 */
    Matrix chol_solve (const Matrix & A, const Matrix & b)
    {
	// NOTE: do not need to check for squareness of A or postive definiteness
	// of A since this is already being done in the cholesky() call below.
	
	if ((b.colsize == 1) && (A.rowsize == b.rowsize) && 
	    (A.rowsize == A.colsize)){
	    
	    register Matrix M = cholesky (A);
	    register double holder;
	    register double *y = new double[A.rowsize];
	    register double *x = new double[A.rowsize];
	    
	    // solve M*y = b
	    for (int i = 0; i < A.rowsize; ++i) {
		holder = 0.0;
		for (int j = 0; j < i; ++j) {
		    holder += M.data[i * M.colsize + j] * y[j];
		}
		y[i] = (1.0 / M.data[i * M.colsize + i]) * (b.data[i] - holder);
	    }
	    
	    // solve M'*x = y
	    for (int i = A.rowsize - 1; i >= 0; --i) {
		holder = 0.0;
		for (int j = i + 1; j < A.rowsize; ++j) {
		    holder += M.data[j * M.colsize + i] * x[j];
		}
		x[i] = (1.0 / M.data[i * M.colsize + i]) * (y[i] - holder);
	    }    
	    
	    Matrix temp = Matrix (x, A.rowsize, 1);
	    delete[]y;
	    delete[]x;
	    return temp;
	} else {
	    cerr <<"ERROR 0021: Inputs not proper dimension in SCYTHE::chol_solve()" 
		 << endl;
	    exit(21);
	}
    }
    
/*! 
 * \brief Solves \a A \a x = \a b for \a x via back-substitution 
 * using Cholesky Decomposition
 *
 * Solves \a A \a x = \a b for \a x via back-substitution using 
 * Cholesky Decomposition.  \a A must be symmetric and positive definite.
 * PLEASE NOTE:  THIS FUNCTION IS OVERLOADED.
 * \remarks The solution is trivial because of the following identity.
 * A*x = (M*M')*x = M*(M'*x) = b.
 * Function works by successively solving two triangular systems:
 * M*y = b for y
 * and then
 * M'*x = y for x.
 \param A a constant reference to the Matrix \a A in the equality
 * \a A \a x = \a b.
 \param b a constant reference to the Matrix \a b in the equality
 * \a A \a x = \a b.
 \param M a constant reference to the Matrix \a M.  See above remarks.
 \return a new Matrix (vector) \e x from the equality
 * \a A \a x = \a b.
 */
    Matrix chol_solve (const Matrix & A, const Matrix & b, const Matrix & M)
    {
	if ((b.colsize == 1) && (A.rowsize == b.rowsize) && 
	    (A.rowsize == M.rowsize) && (A.rowsize == A.colsize) &&
	    (M.rowsize == M.colsize)){
	    
	    register double *y = new double[A.rowsize];
	    register double *x = new double[A.rowsize];
	    
	    // solve M*y = b
	    for (int i = 0; i < A.rowsize; ++i) {
		double holder = 0.0;
		for (int j = 0; j < i; ++j) {
		    holder += M.data[i * M.colsize + j] * y[j];
		}
		y[i] = (1.0 / M.data[i * M.colsize + i]) * (b.data[i] - holder);
	    }
	    
	    // solve M'*x = y
	    for (int i = A.rowsize - 1; i >= 0; --i) {
		double holder = 0.0;
		for (int j = i + 1; j < A.rowsize; ++j) {
		    holder += M.data[j * M.colsize + i] * x[j];
		}
		x[i] = (1.0 / M.data[i * M.colsize + i]) * (y[i] - holder);
	    }
	    
	    Matrix temp = Matrix (x, A.rowsize, 1);
	    delete[]y;
	    delete[]x;
	    return temp;
	} else{
	    cerr << "ERROR 0022: inputs not proper dimension in SCYTHE::chol_solve()" 
		 << endl;
	    exit(22);
	}
    }
    
    
/*! 
 * \brief Calculates the inverse of a Symmetric Positive Definite Matrix 
 *
 * Calculates the inverse of a Symmetric Positive Definite Matrix.  
 * PLEASE NOTE:  THIS FUNCTION IS OVERLOADED.
 * \param A a constant reference to the Matrix to be inverted.
 * \return a new inverted Matrix.
 */
    Matrix invpd (const Matrix & A)
    {
	// SYMMETRY OF A *IS NOT* CHECKED
	register Matrix b = Matrix (A.rowsize, 1);  
	register Matrix Ainv = Matrix(A.rowsize, A.colsize);
	
	// the following block is equivalent to:
	// register Matrix M = cholesky(A);
	if (A.rowsize != A.colsize){
	    cerr << "ERROR 0023: Input Matrix not square in SCYTHE::invpd()\n"
		 << "Cholesky decomposition failed. Exiting program. " 
		 << endl;
	    exit (23);
	}    
	register double *newdata = new double[A.rowsize * A.colsize];
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = i; j < A.colsize; ++j) {
		register double holder = A.data[A.colsize * i + j];
		for (int k = 0; k <= i - 1; ++k) {
		    holder -= newdata[A.colsize * i + k] * newdata[A.colsize * j + k];
		}
		if (i == j) {
		    if (holder <= 0) {
			cerr << 
			    "ERROR 0024: Input Matrix not positive definite in SCYTHE::invpd()\n"
			     << "Cholesky decomposition failed. Exiting program. " 
			     << endl;
			exit (24);
		    }
		    newdata[A.colsize * i + i] = ::sqrt (holder);
		} else{
		    newdata[A.colsize * j + i] = (1.0 / newdata[A.colsize * i + i]) * holder;
		    newdata[A.colsize * i + j] = 0.0;
		}
	    }
	}    
	register Matrix M = Matrix (newdata, A.rowsize, A.colsize);
	delete [] newdata;
	// end register Matrix M = cholesky(A); block 
	
	// following 2 lines are actually part of the 
	// cholesky solve block below
	register double *y = new double[A.rowsize];
	register double *x = new double[A.rowsize];
	
	for (int j = 0; j < A.rowsize; ++j) {
	    b.data[j] = 1.0;
	    // The following block is equivalent to: 
	    // Matrix hold = chol_solve (A, b, M);
	    
	    // solve M*y = b
	    for (int i = 0; i < A.rowsize; ++i) {
		double holder = 0.0;
		for (int k = 0; k < i; ++k) {
		    holder += M.data[i * M.colsize + k] * y[k];
		}
		y[i] = (1.0 / M.data[i * M.colsize + i]) * (b.data[i] - holder);
	    }
	    
	    // solve M'*x = y 
	    for (int i = A.rowsize - 1; i >= 0; --i) {
		double holder = 0.0;
		for (int k = i + 1; k < A.rowsize; ++k) {
		    holder += M.data[k * M.colsize + i] * x[k];
		}
		x[i] = (1.0 / M.data[i * M.colsize + i]) * (y[i] - holder);
	    }
	    Matrix hold = Matrix (x, A.rowsize, 1);      
	    // end Matrix hold = chol_solve(A, b, M); block
	    
	    b.data[j] = 0.0;
	    for (int k=0; k<A.rowsize; ++k)
		Ainv(k,j) = hold[k];
	}
	
	delete[]y;
	delete[]x;
	
	return Ainv;
    }
    
/*!
 * \brief Calculates the inverse of a Symmetric Positive Definite Matrix 
 *
 * Calculates the inverse of a Symmetric Positive Definite Matrix.  
 * PLEASE NOTE:  THIS FUNCTION IS OVERLOADED.  
 * SYMMETRY OF \a A IS NOT CHECKED.
 * \param A a constant reference to the Matrix to be inverted.
 * \param M a constant reference to a Matrix from the Cholesky 
 * Decomposition of Matrix /a A.
 * \return a new inverted Matrix.
 */
    Matrix invpd (const Matrix & A, const Matrix & M)
    {
	
	register Matrix b = Matrix (A.rowsize, 1); 
	register Matrix Ainv = Matrix(A.rowsize, A.rowsize);    
	
	// following 2 lines are actually part of the 
	// cholesky solve block below
	register double *y = new double[A.rowsize];
	register double *x = new double[A.rowsize];
	
	for (int j = 0; j < A.rowsize; ++j) {
	    b.data[j] = 1.0;
	    
	    // The following block is equivalent to: 
	    // Matrix hold = chol_solve (A, b, M);
	    
	    // solve M*y = b
	    for (int i = 0; i < A.rowsize; ++i) {
		double holder = 0.0;
		for (int k = 0; k < i; ++k) {
		    holder += M.data[i * M.colsize + k] * y[k];
		}
		y[i] = (1.0 / M.data[i * M.colsize + i]) * (b.data[i] - holder);
	    }
	    
	    // solve M'*x = y 
	    for (int i = A.rowsize - 1; i >= 0; --i) {
		double holder = 0.0;
		for (int k = i + 1; k < A.rowsize; ++k) {
		    holder += M.data[k * M.colsize + i] * x[k];
		}
		x[i] = (1.0 / M.data[i * M.colsize + i]) * (y[i] - holder);
	    }
	    Matrix hold = Matrix (x, A.rowsize, 1);      
	    // end Matrix hold = chol_solve(A, b, M); block
	    
	    b.data[j] = 0.0;
	    for (int k=0; k<A.rowsize; ++k)
		Ainv(k,j) = hold[k];
	}
	
	delete[]y;
	delete[]x;
	
	return Ainv;
    }
    
//  This code is based on  Algorithm 3.4.1 of Golub and Van Loan 
//  3rd edition, 1996. Major difference is in how the output is 
//  structured. 
    
/*! 
 * \brief Calculates the LU Decomposition of a square Matrix 
 *
 * Calculates the LU Decomposition of a square Matrix.
 * \remarks This function solves \a P\a A \a = \a L \a U
 * for \a L and \a U where \a P is a permutation Matrix and \a A 
 * is the input Matrix also returns the permutation vector as 
 * defined in Golub and Van Loan.
 * \param AA a constant reference to the Matrix \a A (see above).
 * \param L a reference to the Matrix \a L (see above).
 * \param U a reference to the Matrix \a U (see above).
 * \param perm_vec a references to the Matrix \a P (see above).
 * \return an integer ( 0 if successful, otherwise it returns 
 * an error number).
 */
    int lu_decomp(const Matrix& AA, Matrix& L, Matrix& U, Matrix& perm_vec){
	Matrix A = AA;
	
	if (A.rowsize != A.colsize) {
	    cerr << "ERROR 0025: Matrix A not square in SCYTHE::lu_decomp()" 
		 << endl;
	    exit(25);
	}
	
	if (A.rowsize == 1) {
	    L = ones(1,1);
	    U = AA;
	    perm_vec = Matrix(1,1);
	    return 0;
	}
	
	L = U = Matrix(A.rowsize, A.rowsize);
	perm_vec = Matrix(A.rowsize-1,1);
	
	for (int k=0; k<A.rowsize-1; ++k) {
	    int pivot = k;
	    // find pivot
	    for (int i=k; i<A.rowsize; ++i){
		if ( ::fabs(A(pivot,k)) < ::fabs(A(i,k)) ) {
		    pivot = i; 
		} 
	    }
	    
	    if(A(pivot,k) == 0.0) {
		cerr << "ERROR 0026: Matrix A is singular in SCYTHE::lu_decomp()" 
		     << endl;
		exit(26);
	    }
	    
	    // permute 
	    if (k != pivot){
		for (int i=0; i<A.rowsize; ++i){
		    double temp = A(pivot,i);
		    A(pivot,i) = A(k,i);
		    A(k,i) = temp;
		}
	    }
	    perm_vec[k] = pivot;
	    
	    for (int i = k+1; i<A.rowsize; ++i){
		A(i,k) = A(i,k)/A(k,k);
		for (int j = k+1; j <A.rowsize; ++j){
		    A(i,j) = A(i,j) - A(i,k)*A(k,j);
		}
	    }
	}
	
	L = A; 
	for (int i=0; i<A.rowsize; ++i){
	    for (int j=i; j<A.rowsize; ++j){
		U(i,j) = A(i,j);
		L(i,j) = 0.0;
		L(i,i) = 1.0;
	    }
	}
	
	return 0;
    }
    
    
    
//  4/29/2001 (KQ)
/*! 
 * \brief Solves \a A*x=b for \a x
 *
 * Solves \a A \a x \a = \a b for \a x via forward and 
 * back-substitution using the LU Decomposition of \a A.  
 * PLEASE NOTE: THIS FUNCTION IS OVERLOADED.
 * \param AA a constant reference to the Matrix \a A.
 * \param b a reference to the Matrix \a b.
 * \return a Matrix object \a x.
 */
    Matrix lu_solve(const Matrix& AA, const Matrix& b){
	if (b.colsize != 1){
	    cerr << "Error 0027: Vector b not column vector in SCYTHE::lu_solve()" 
		 << endl;
	    exit(27);
	}
	if (AA.rowsize != AA.colsize){
	    cerr << "Error 0028: Matrix A not square in SCYTHE::lu_solve()" 
		 << endl;
	    exit(28);
	}
	if (AA.rowsize != b.rowsize){
	    cerr << "Error 0029: Matrix A and b not conformable in SCYTHE::lu_solve()" 
		 << endl;
	    exit(29);
	}
	
	// step 1 compute the LU factorization 
	// taken from lu_decomp()
	Matrix A = AA;
	Matrix L, U, perm_vec;
	
	if (A.rowsize == 1){
	    L = ones(1,1);
	    U = AA;
	    perm_vec = Matrix(1,1);
	} else {
	    L = U = Matrix(A.rowsize, A.rowsize);
	    perm_vec = Matrix(A.rowsize-1,1);
	    
	    for (int k=0; k<A.rowsize-1; ++k){
		int pivot = k;
		// find pivot
		for (int i=k; i<A.rowsize; ++i){
		    if ( ::fabs(A(pivot,k)) < ::fabs(A(i,k)) ) pivot = i;  
		}
		
		if(A(pivot,k) == 0.0){
		    cerr << "ERROR 0030: Matrix A is singular in SCYTHE::lu_solve()" 
			 << endl;
		    exit(30);
		}
		
		// permute 
		if (k != pivot){
		    for (int i=0; i<A.rowsize; ++i){
			double temp = A(pivot,i);
			A(pivot,i) = A(k,i);
			A(k,i) = temp;
		    }
		}
		perm_vec[k] = pivot;
		
		for (int i = k+1; i<A.rowsize; ++i){
		    A(i,k) = A(i,k)/A(k,k);
		    for (int j = k+1; j <A.rowsize; ++j){
			A(i,j) = A(i,j) - A(i,k)*A(k,j);
		    }
		}
	    }
	    
	    L = A; 
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=i; j<A.rowsize; ++j){
		    U(i,j) = A(i,j);
		    L(i,j) = 0.0;
		    L(i,i) = 1.0;
		}
	    }
	}
	
	// step 2 solve L*y = Pb via forward substitution
	Matrix bb = row_interchange(b, perm_vec);
	Matrix y = Matrix(A.rowsize,1);
	for (int i=0; i<A.rowsize; ++i){
	    double sum = 0.0;
	    for (int j=0; j<i; ++j)
		sum += L.data[j+A.colsize*i]*y.data[j];
	    y.data[i] = (bb.data[i] - sum)/L.data[i+A.colsize*i];
	}
	
	// step 3 solve U*x = y via backsubstitution
	Matrix x = Matrix(A.rowsize,1);
	for (int i=A.rowsize-1; i>=0; --i){
	    double sum = 0.0;
	    for (int j=i+1; j<A.rowsize; ++j)
		sum += U.data[j+A.colsize*i] * x.data[j];
	    x.data[i] = (y.data[i] - sum)/U.data[i+A.colsize*i];
	}
	
	return x;
    }
    
//  4/29/2001 (KQ)
/*! 
 * \brief Solves \a A*x=b for \a x
 *
 * Solves \a A \a x \a = \a b for \a x via forward and 
 * back-substitution using the LU Decomposition of \a A.  
 * PLEASE NOTE: THIS FUNCTION IS OVERLOADED.
 * \param A a constant reference to the Matrix \a A.
 * \param b a constant reference to the Matrix \a b.
 * \param L a constant reference to the Matrix \a L used in the LU 
 * Decomposition.
 * \param b a constant reference to the Matrix \a U used in the LU 
 * Decomposition.
 * \param p a vector used in the LU Decomposition.
 * \return a Matrix object \a x.
 */
    Matrix lu_solve(const Matrix& A, const Matrix& b, 
		    const Matrix& L, const Matrix& U, const Matrix& p){
	if (b.colsize != 1){
	    cerr << "Error 0031: Vector b not column vector in SCYTHE::lu_solve()" 
		 << endl;
	    exit(31);
	}
	if (A.rowsize != A.colsize){
	    cerr << "Error 0032: Matrix A not square in SCYTHE::lu_solve()" 
		 << endl;
	    exit(32);
	}
	if (A.rowsize != b.rowsize){
	    cerr << "Error 0033: Matrices A and b not conformable in SCYTHE::lu_solve()" 
		 << endl;
	    exit(33);
	}
	if ((A.rowsize != L.rowsize) || (A.rowsize != U.rowsize) ||
	    (A.colsize != L.colsize) || (A.colsize != U.colsize)){
	    cerr << "Error 0034: Matrices A, L, and U not of same dimension in "
		 << "SCYTHE::lu_solve()" << endl;
	    exit(34);
	}
	if ( (p.rowsize +1) != A.rowsize){
	    cerr << "ERROR 0035: Matrices A and p not of consistent sizes in "
		 << "SCYTHE::lu_solve()"  << endl;
	    exit(35);
	}
	
	// step 1 solve L*y = Pb via forward substitution
	Matrix bb = row_interchange(b, p);
	Matrix y = Matrix(A.rowsize,1);
	for (int i=0; i<A.rowsize; ++i){
	    double sum = 0.0;
	    for (int j=0; j<i; ++j)
		sum += L.data[j+A.colsize*i]*y.data[j];
	    y.data[i] = (bb.data[i] - sum)/L.data[i+A.colsize*i];
	}
	
	// step 2 solve U*x = y via backsubstitution
	Matrix x = Matrix(A.rowsize,1);
	for (int i=A.rowsize-1; i>=0; --i){
	    double sum = 0.0;
	    for (int j=i+1; j<A.rowsize; ++j)
		sum += U.data[j+A.colsize*i] * x.data[j];
	    x.data[i] = (y.data[i] - sum)/U.data[i+A.colsize*i];
	}
	
	return x;
    }
    
    
//  4/29/2001 (KQ)
/*! 
 * \brief Interchanges the rows of \a A with those in vector \a p
 *
 * Interchanges the rows of \a A with those in vector \a p and 
 * returns the modified Matrix. Useful for putting \a A into the form 
 * of its permuted LU factorization.
 * \param A a constant reference to the Matrix \a A.
 * \param pp a constant reference to the Matrix (vector) \a p from 
 * which the interchange row will come from.
 * \return the modified Matrix \a A.
 */
    Matrix row_interchange(const Matrix& A, const Matrix& pp){
	Matrix PA = A;
	Matrix p = pp;
	if (p.colsize != 1){
	    cerr << "ERROR 0036: Vector p not a column vector in "
		 << "SCYTHE::row_interchange()" << endl;
	    exit(36);
	}
	if ( (p.rowsize +1) != A.rowsize){
	    cerr << "ERROR 0037: Matrices A and p not of consistent sizes in "
		 << "SCYTHE::row_interchange()"  << endl;
	    exit(37);
	}
	
	for (int i=0; i<(A.rowsize-1); ++i){
	    //swap A(i,.) and A(p[i],.)
	    int swap_row = static_cast<int>(p.data[i]);
	    for (int j=0; j<A.colsize; ++j){
		double temp = PA.data[j+A.colsize*i];
		PA.data[j+A.colsize*i] = PA.data[j+A.colsize*swap_row];
		PA.data[j+A.colsize*swap_row] = temp;
	    }
	}
	
	return PA;
    }
    
    
//  4/29/2001 (KQ)
/*! 
 * \brief Calculate the Inverse of a square Matrix \a A
 *
 * Calculate the Inverse of a square Matrix \a A via LU decomposition.
 * \param AA a constant reference to the Matrix \a A to be inverted.
 * \return the inverted Matrix \a A.
 */
    Matrix inv (const Matrix & AA)
    {
	if (AA.rowsize != AA.colsize){
	    cerr << "Error 0038: Matrix A not square in SCYTHE::inv()" 
		 << endl;
	    exit(38);
	}
	
	Matrix b = Matrix (AA.rowsize, 1); 
	Matrix Ainv = Matrix(AA.rowsize, AA.rowsize);    
	
	// step 1 compute the LU factorization 
	// taken from lu_decomp()
	Matrix A = AA;
	Matrix L, U, perm_vec;
	
	if (A.rowsize == 1){
	    L = ones(1,1);
	    U = AA;
	    perm_vec = Matrix(1,1);
	} else {
	    L = U = Matrix(A.rowsize, A.rowsize);
	    perm_vec = Matrix(A.rowsize-1,1);
	    
	    for (int k=0; k<A.rowsize-1; ++k){
		int pivot = k;
		// find pivot
		for (int i=k; i<A.rowsize; ++i){
		    if ( ::fabs(A(pivot,k)) < ::fabs(A(i,k)) ) pivot = i;  
		}
		
		if(A(pivot,k) == 0.0){
		    cerr << "ERROR 0039: Matrix A is singular in SCYTHE::inv()" 
			 << endl;
		    exit(39);
		}
		
		// permute 
		if (k != pivot){
		    for (int i=0; i<A.rowsize; ++i){
			double temp = A(pivot,i);
			A(pivot,i) = A(k,i);
			A(k,i) = temp;
		    }
		}
		perm_vec[k] = pivot;
		
		for (int i = k+1; i<A.rowsize; ++i){
		    A(i,k) = A(i,k)/A(k,k);
		    for (int j = k+1; j <A.rowsize; ++j){
			A(i,j) = A(i,j) - A(i,k)*A(k,j);
		    }
		}
	    }
	    
	    L = A; 
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=i; j<A.rowsize; ++j){
		    U(i,j) = A(i,j);
		    L(i,j) = 0.0;
		    L(i,i) = 1.0;
		}
	    }
	}
	
	// step 2 repeated solving of A*hold = b
	for (int j = 0; j < A.rowsize; ++j) {
	    b.data[j] = 1.0;
	    //Matrix hold = lu_solve(A, b, L, U, p);
	    
	    // step 2.1 solve L*y = Pb via forward substitution
	    Matrix bb = row_interchange(b, perm_vec);
	    Matrix y = Matrix(A.rowsize,1);
	    for (int i=0; i<A.rowsize; ++i){
		double sum = 0.0;
		for (int j=0; j<i; ++j)
		    sum += L.data[j+A.colsize*i]*y.data[j];
		y.data[i] = (bb.data[i] - sum)/L.data[i+A.colsize*i];
	    }
	    
	    // step 2.2 solve U*x = y via backsubstitution
	    Matrix x = Matrix(A.rowsize,1);
	    for (int i=A.rowsize-1; i>=0; --i){
		double sum = 0.0;
		for (int j=i+1; j<A.rowsize; ++j)
		    sum += U.data[j+A.colsize*i] * x.data[j];
		x.data[i] = (y.data[i] - sum)/U.data[i+A.colsize*i];
	    }
	    
	    // step 3 reset b and put the solution in Ainv
	    b.data[j] = 0.0;
	    for (int k=0; k<A.rowsize; ++k)
		Ainv(k,j) = x[k];
	}
	
	return Ainv;
    }
    
    
//  NOTE: LU decomposition algorithm is based on  Algorithm 3.4.1 
//      of Golub and Van Loan 3rd edition, 1996. 
/*! 
 * \brief Calculate the determinant of a square Matrix \a A
 *
 * Calculate the determinant of a square Matrix \a A via LU decomposition.
 * \param AA a constant reference to the Matrix.
 * \return the determinant of \a A.
 */
    double det(const Matrix& AA){
	
	Matrix A = AA;
	
	if (A.rowsize != A.colsize){
	    cerr << "ERROR 0040: Matrix A not square in SCYTHE::det()" << endl;
	    exit(40);
	}
	
	if(A.rowsize == 1)
	    return A(0,0);
	
	Matrix L = Matrix(A.rowsize, A.rowsize);
	Matrix U = L;
	double sign = 1.0;
	
	for (int k=0; k<A.rowsize-1; ++k){
	    int pivot = k;
	    // find pivot
	    for (int i=k; i<A.rowsize; ++i){
		if ( A(pivot,k) < ::fabs(A(i,k)) ) pivot = i;  
	    }
	    
	    if(A(pivot,k) == 0.0){
		return 0.0;
	    }
	    
	    // permute 
	    if (k != pivot){
		sign = -1*sign;
		for (int i=k; i<A.rowsize; ++i){
		    double temp = A(pivot,i);
		    A(pivot,i) = A(k,i);
		    A(k,i) = temp;
		}
	    }
	    
	    for (int i = k+1; i<A.rowsize; ++i){
		A(i,k) = A(i,k)/A(k,k);
		for (int j = k+1; j <A.rowsize; ++j){
		    A(i,j) = A(i,j) - A(i,k)*A(k,j);
		}
	    }
	}
	
	double det = 1.0;
	for (int i = 0; i<A.rowsize; ++i)
	    det = det*A(i,i);
	
	return sign*det;
    }
    
//! Column bind 2 matrices
/*!
 * Column bind 2 matrices,\a A and \a B.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the new Matrix \a A
 */
    Matrix cbind (const Matrix & A, const Matrix & B)
    {
	if (A.rowsize != B.rowsize) {
	    cerr << "ERROR 0041: Matrices A and B do not have some number of "
		 << "rows in SCYTHE::cbind(). \nSCYTHE::cbind() failed. Exiting program." 
		 << endl;
	    exit (41);
	}
	
	int totalcols = A.colsize + B.colsize;
	double *newdata = new double[A.rowsize * totalcols];
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[i * totalcols + j] = A.data[i * A.colsize + j];
	    }
	    for (int k = 0; k < B.colsize; ++k) {
		newdata[i * totalcols + k + A.colsize] = B.data[i * B.colsize + k];
	    }
	}
	
	Matrix temp = Matrix (newdata, A.rowsize, totalcols);
	delete[]newdata;
	return temp;
    }
    
//! Row bind 2 matrices
/*!
 * Row bind 2 matrices,\a A and \a B.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the new Matrix \a A.
 */
    Matrix rbind (const Matrix & A, const Matrix & B)
    {
	if (A.colsize != B.colsize) {
	    cerr << "ERROR 0042: Matrices A and B do not have some number of "
		 << "cols in SCYTHE::rbind(). \nSCYTHE::rbind() failed. Exiting program." 
		 << endl;
	    exit (42);
	}
	
	int totalrows = A.rowsize + B.rowsize;
	double *newdata = new double[totalrows * A.colsize];
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j)  {
		newdata[i * A.colsize + j] = A.data[i * A.colsize + j];
	    }
	}
	for (int k = 0; k < B.rowsize; ++k) {
	    for (int j = 0; j < B.colsize; ++j) {
		newdata[k * B.colsize + (A.rowsize * A.colsize) + j] =
		    B.data[k * B.colsize + j];
	    }
	}
	
	Matrix temp = Matrix (newdata, totalrows, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the sum of each column of a Matrix
/*!
 * Calculate the sum of each column of a Matrix.  
 * Returns the row vector of sums.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of sums for each corresponding column.
 */
    Matrix sumc (const Matrix & A)
    {
	double *newdata = new double[A.colsize];
	
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[j] += A.data[A.colsize * i + j];
	    }
	}
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the product of each column of a Matrix
/*!
 * Calculate the product of each column of a Matrix.  
 * Returns the row vector of products.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of products for each corresponding column.
 */
    Matrix prodc (const Matrix & A)
    {
	double *newdata = new double[A.colsize];
	
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = 1.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[j] = newdata[j] * A.data[A.colsize * i + j];
	    }
	}
	
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the mean of each column of a Matrix
/*!
 * Calculate the mean of each column of a Matrix.  
 * Returns the row vector of means.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of means for each corresponding column.
 */
    Matrix meanc (const Matrix & A)
    {
	double *newdata = new double[A.colsize];
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[j] += A.data[A.colsize * i + j];
	    }
	}
	for (int i = 0; i < A.colsize; ++i) {
	    newdata[i] = (1.0 / A.rowsize) * newdata[i];
	}
	
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the variances of each column of a Matrix
/*!
 * Calculate the variances of each column of a Matrix.  
 * Returns the row vector of variances.
 * \param A a constant reference to a Matrix \a A.
 * \return the vector of variances for each corresponding column.
 */
    Matrix varc (const Matrix & A)
    {
	Matrix mu = meanc (A);
	double *newdata = new double[A.colsize];
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[j] += ::pow (mu.data[j] - A.data[A.colsize * i + j], 2) *
		    (1.0 / (A.rowsize-1));
	    }
	}
	
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the standard deviation of each column of a Matrix
/*!
 * Calculate the standard deviation of each column of a Matrix.  
 * Returns the row vector of standard deviations.
 * \param A a constant reference to a Matrix \aA.
 * \return the vector of standard deviations for each corresponding column.
 */
    Matrix stdc (const Matrix & A)
    {
	Matrix mu = meanc (A);
	double *newdata = new double[A.colsize];
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = 0.0;
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		newdata[j] += ::pow (mu.data[j] - A.data[A.colsize * i + j], 2) *
		    (1.0 / (A.rowsize-1));
	    }
	}
	
	for (int i = 0; i < A.colsize; ++i)
	    newdata[i] = ::sqrt (newdata[i]);
	
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the square root of each element of a Matrix
/*!
 * Calculate the square root of each element of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the Matrix of square roots.
 */
    Matrix sqrt (const Matrix & A)
    {
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::sqrt (A.data[i]);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the absolute value of each element of a Matrix
/*!
 * Calculate the absolute value of each element of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the Matrix of absolute values.
 */
    Matrix fabs (const Matrix & A)
    {
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::fabs (A.data[i]);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
    
//! Calculate the value of \a e^x for each individual Matrix element
/*!
 * Calculate the value of \a e^x for each individual Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of exponentials.
 */
    Matrix exp(const Matrix& A){
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::exp (A.data[i]);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the natural log of each individual Matrix element
/*!
 * Calculate the natural log of each individual Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of natural logs.
 */
    Matrix log(const Matrix& A){
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::log (A.data[i]);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculate the Base 10 Log of each Matrix element
/*!
 * Calculate the Base 10 Log of each Matrix element.
 * \param A a constant reference to a Matrix \a A.
 * \return the new Matrix of Base 10 logs.
 */
    Matrix log10(const Matrix& A){
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::log10(A.data[i]);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Raises each Matrix element to a specified power
/*!
 * Raises each Matrix element to a specified power.
 * \param A a constant reference to a Matrix \a A.
 * \return the new modified Matrix.
 */
    Matrix pow(const Matrix& A, const double& e){
	double *newdata = new double[A.size];
	
	for (int i = 0; i < A.size; ++i)
	    newdata[i] = ::pow(A.data[i], e);
	
	Matrix temp = Matrix (newdata, A.rowsize, A.colsize);
	delete[]newdata;
	return temp;
    }
    
    
//! Calculates the maximum element in a Matrix
/*!
 * Calculates the maximum element in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the maximum element (a double).
 */
    double max (const Matrix & A)
    {
	double max = A.data[0];
	for (int i = 1; i < A.size; ++i) {
	    if (A.data[i] > max)
		max = A.data[i];
	}
	return max;
    }
    
//! Calculates the minimum element in a Matrix
/*!
 * Calculates the minimum element in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return the minimum element (a double).
 */
    double min (const Matrix & A)
    {
	double min = A.data[0];
	for (int i = 1; i < A.size; ++i) {
	    if (A.data[i] < min)
		min = A.data[i];
	}
	return min;
    }
    
//! Calculates the maximum of each Matrix column
/*!
 * Calculates the maximum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the maximum elements.
 */
    Matrix maxc (const Matrix & A)
    {
	double *newdata = new double[A.colsize];
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		if (i == 0) {
		    newdata[j] = A.data[A.colsize * i + j];
		} else if (A.data[A.colsize * i + j] > newdata[j]) {
		    newdata[j] = A.data[A.colsize * i + j];
		}
	    }
	}
	
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
//! Calculates the minimum of each Matrix column
/*!
 * Calculates the minimum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the minimum elements.
 */
    Matrix minc (const Matrix & A)
    {
	double *newdata = new double[A.colsize];
	
	for (int i = 0; i < A.rowsize; ++i) {
	    for (int j = 0; j < A.colsize; ++j) {
		if (i == 0) {
		    newdata[j] = A.data[A.colsize * i + j];
		} else if (A.data[A.colsize * i + j] < newdata[j]) {
		    newdata[j] = A.data[A.colsize * i + j];
		}
	    }
	}
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
    
//! Finds the index of the maximum of each Matrix column
/*!
 * Finds the index of the maximum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the index of each maximum element.
 */
    Matrix maxindc(const Matrix& A){
	double *newdata = new double[A.colsize];
	Matrix maxvec = Matrix(1,A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i){
	    for (int j = 0; j < A.colsize; ++j){
		if (i == 0){
		    maxvec[j] = A.data[A.colsize * i + j]; 
		    newdata[j] = 0;
		} else if (A.data[A.colsize * i + j] > maxvec[j]){
		    maxvec[j] = A.data[A.colsize * i + j];
		    newdata[j] = i;
		}
	    }
	}
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
    
//! Finds the index of the minimum of each Matrix column
/*!
 * Finds the index of the minimum of each Matrix column.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) of the index of each minimum element.
 */
    Matrix minindc(const Matrix& A){
	double *newdata = new double[A.colsize];
	Matrix minvec = Matrix(1,A.colsize);
	
	for (int i = 0; i < A.rowsize; ++i){
	    for (int j = 0; j < A.colsize; ++j){
		if (i == 0){
		    minvec[j] = A.data[A.colsize * i + j]; 
		    newdata[j] = 0;
		} else if (A.data[A.colsize * i + j] < minvec[j]){
		    minvec[j] = A.data[A.colsize * i + j];
		    newdata[j] = i;
		}
	    }
	}
	Matrix temp = Matrix (newdata, 1, A.colsize);
	delete[]newdata;
	return temp;
    }
    
    
//! Calculates the order of each element in a Matrix
/*!
 * Calculates the order of each element in a Matrix.
 * \param A a constant reference to a Matrix A.
 * \return a Matrix (vector) in which the \e i'th element 
 * gives the order position of the \e i'th element of \a A.
 */
    Matrix order(const Matrix& A){
	if (A.colsize != 1){
	    cerr << "ERROR 0043: Matrix A not a column vector in SCYTHE::order()" 
		 << endl;
	    exit(43);
	}
	double *newdata = new double[A.rowsize];
	
	for (int i=0; i<A.rowsize; ++i){
	    newdata[i] = sumc(A << A.data[i])[0];    
	}
	
	Matrix temp = Matrix (newdata, A.rowsize, 1);
	delete[] newdata;
	return temp;
    }
    
    
//! Selects all the rows of Matrix \a A for which binary column vector \a e has an element equal to 1
/*!
 * Selects all the rows of Matrix \a A for which binary column 
 * vector \a e has an element equal to 1.
 * \param A a constant reference to a Matrix \a A (n x k).
 * \param e a constant reference to a vector \a e (n x 1).
 * \return a Matrix of all rows of \a A for \a e equal to one.
 */
    Matrix selif(const Matrix& A, const Matrix& e){
	
	// check to see if rowsize matches  
	if (A.rowsize != e.rowsize){
	    cerr << "ERROR 0044: Matrices not conformable in SCYTHE::selif()" 
		 << endl;
	    exit(44);
	}
	
	// check to see if e is a column vector
	if (e.colsize > 1){
	    cerr << "ERROR 0045: e not a column vector in SCYTHE::selif()" 
		 << endl;
	    exit(45);
	}
	
	// loop to check if e contains binary data, and count number
	// of output rows 
	int N = 0;
	for (int i=0; i<e.rowsize; ++i){
	    if (e.data[i] != 0 && e.data[i] != 1){
		cerr << "ERROR 0046: Vector e contains non binary data in "
		     <<" SCYTHE::selif()" 
		     << endl;
		exit(46);
	    }
	    if (e.data[i] == 1) {
		N += 1;
	    }
	}
	
	// declare and form output Matrix
	double *newdata = new double[A.colsize*N];
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    if (e.data[i] == 1){
		for (int j=0; j<A.colsize; ++j){
		    newdata[count] = A.data[A.colsize*i + j]; 
		    ++count;
		}
	    }
	}
	Matrix temp = Matrix(newdata, N, A.colsize);
	delete[] newdata;
	return temp;
    }
    
    
//! Finds unique elements in a Matrix
/*!
 * Finds unique elements in a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix of all unique elements.
 */
    Matrix unique(const Matrix& A){
	double *newdata = new double[A.size];
	
	newdata[0] = A.data[0];
	int count = 1;
	for (int i=1; i<A.size; ++i){
	    int uniq = 1;
	    for (int j=0; j<count; ++j){
		if (newdata[j] == A.data[i]){
		    uniq = 0;
		    break;
		}
	    }
	    if (uniq==1){
		newdata[count] = A.data[i];
		++count;
	    }
	}
	
	Matrix temp = Matrix(newdata, count, 1);
	delete[] newdata;
	return temp;
    }
    
//! Turn Matrix into Column vector by stacking columns
/*!
 * Turn Matrix into Column vector by stacking columns. 
 * NOTE: \e Vecr() is much faster than \e Vecc().
 * \param A a constant reference to a Matrix \a A.
 * \return a Column vector.
 */
    Matrix vecc(const Matrix& A){
	// first transposes the input Matrix's data 
	int newrowsize = A.colsize;
	int newcolsize = A.rowsize;
	double *newdata = new double[A.size];
	for (int i = 0; i < newcolsize; ++i) {
	    for (int j = 0; j < newrowsize; ++j) {
		newdata[i + newcolsize * j] = A.data[j + newrowsize * i];
	    }
	}
	// then takes this transposed data and vecrs
	Matrix temp = Matrix (newdata, A.size, 1);
	delete[]newdata;
	return temp;
    }
    
    
//! Reshapes a row major order Matrix or Vector
/*!
 * Reshapes a row major order Matrix or Vector.
 * \param A a constant reference to a Matrix \a A.
 * \param r a constant integer reflecting the new number of rows.
 * \param c a constant integer reflecting the new number of columns.
 * \return a new Matrix with the specified dimensions.
 */
    Matrix reshape(const Matrix& A, const int r, const int c){
	if (A.size != r*c){
	    cerr << "ERROR 0047: Input dimensions to SCYTHE::reshape() not consistent "
		 << "with size of input Matrix" 
		 << endl;
	    exit(47);
	}
	Matrix temp = Matrix(A.data, r, c);
	return temp;
    }
    
//! Make vector out of unique elements of a symmetric Matrix
/*!
 * Make vector out of unique elements of a symmetric Matrix.  
 * NOTE: DOES NOT CHECK FOR SYMMETRY!!!
 * \param A a constant reference to a Matrix \a A.
 * \return a Column vector.
 */
    Matrix vech(const Matrix& A){
	if (A.rowsize != A.colsize){
	    cerr << "ERROR 0048: Input Matrix not square in SCYTHE::vech()" << endl;
	    exit(48);
	}
	
	int newsize = static_cast<int>(0.5*(A.size - A.rowsize) + A.rowsize);
	double* newdata = new double[newsize];
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    for(int j=i; j<A.colsize; ++j){
		newdata[count] = A.data[i*A.colsize + j];
		++count;
	    }
	}
	
	Matrix temp = Matrix(newdata, newsize, 1);
	delete [] newdata;
	return temp;
    }
    
    
//! Get symmetric Matrix B back from \a A = \a vech(B)
/*!
 * Get symmetric Matrix \a B back from \a A = \a vech(B)
 * \param A a constant reference to a Matrix \a A.
 * \return a Symmetric Matrix.
 */
    Matrix xpnd(const Matrix& A){
	double newrowsize_d = -.5 + .5*::sqrt(1+8*A.size);
	if (fmod(newrowsize_d,1.0) != 0.0){
	    cerr << "ERROR 0049: Not possible to make square Matrix out of "
		 << "input Matrix to SCYTHE::xpnd()" 
		 << endl;
	    exit(49);
	}
	int newrowsize = static_cast<int>(newrowsize_d);
	double* newdata = new double[newrowsize*newrowsize];
	int count = 0;
	for (int i=0; i<newrowsize; ++i){
	    for (int j=i; j<newrowsize; ++j){
		newdata[i*newrowsize +j] = newdata[j*newrowsize + i] = A.data[count];
		++count;
	    }
	}
	Matrix temp = Matrix(newdata, newrowsize, newrowsize);
	delete [] newdata;
	return temp;
    }
    
    
//! Get the diagonal of a Matrix
/*!
 * Get the diagonal of a Matrix.
 * \param A a constant reference to a Matrix \a A.
 * \return a Matrix (vector) containing the diagonal of \a A.
 */
    Matrix diag(const Matrix& A){
	if (A.rowsize != A.colsize){
	    cerr << "ERROR 0050: Matrix A not square in SCYTHE::diag()" 
		 << endl;
	    exit(50);
	}
	double* newdata = new double[A.rowsize];
	for (int i=0; i<A.rowsize; ++i){
	    newdata[i] = A.data[i*A.colsize + i];
	}
	Matrix temp = Matrix(newdata, A.rowsize, 1);
	delete [] newdata;
	return temp;
    }
    
    
//! Fast calculation of \a A * \a B + \a C
/*!
 * Fast calculation of \a A * \a B + \a C
 * \param A a constant reference to a Matrix \a A.
 * \return the final Matrix (the value of \a A * \a B + \a C).
 */
    Matrix gaxpy(const Matrix& A, const Matrix& B, const Matrix& C){
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	    if (B.rowsize == C.rowsize && B.colsize == C.colsize){
		double *prod = new double[B.size];
		for (int i = 0; i < B.size; ++i) {
		    prod[i] = A.data[0] * B.data[i] + C.data[i];
		}
		Matrix temp = Matrix (prod, B.rowsize, B.colsize);
		delete[]prod;
		return temp;
	    } else {
		cerr << "ERROR 0051: A*B and C not conformable in SCYTHE::gaxpy()" 
		     << endl;
		exit(51);
	    }
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    if (A.rowsize == C.rowsize && A.colsize == C.colsize){
		double *prod = new double[A.size];
		for (int i = 0; i < A.size; ++i) {
		    prod[i] = A.data[i] * B.data[0] + C.data[i];
		}
		Matrix temp = Matrix (prod, A.rowsize, A.colsize);
		delete[]prod;
		return temp;
	    } else {
		cerr << "ERROR 0052: A*B and C not conformable in SCYTHE::gaxpy()" 
		     << endl;
		exit(52);
	    }
	} else if (A.colsize != B.rowsize) {
	    // Case 3: A is n x k and B is m x j (m !=j)
	    cerr << "ERROR 0053: Matrices not conformable for multiplication "
		 << "in SCYTHE::gaxpy()\n exiting program due to error" 
		 << endl;
	    exit (53);
	} else if (A.rowsize == C.rowsize && B.colsize == C.colsize){
	    // Case 4: A is n x k and B is k x j
	    register double *newdata = new double[A.rowsize * B.colsize];
	    for (int i = 0; i < A.rowsize; ++i){
		for (int j = 0; j < B.colsize; ++j){
		    newdata[i*B.colsize + j] = C.data[i*B.colsize + j];
		    for (int k = 0; k < B.rowsize; ++k){
			newdata[i * B.colsize + j] += A.data[i * A.colsize + k] *
			    B.data[k * B.colsize + j];
		    }
		}
	    }
	    
	    Matrix temp = Matrix (newdata, A.rowsize, B.colsize);
	    delete[]newdata;
	    return temp;
	} else {
	    cerr << "ERROR 0054: A*B and C not conformable in SCYTHE::gaxpy()" 
		 << endl;
	    exit(54);
	}
    }
    
    
//! Fast calculation of \a A'A
/*!
 * Fast calculation of \a A'A
 * \param A a constant reference to a Matrix \a A.
 * \return the final Matrix (the value of \a A'A). 
 */
// original
    Matrix crossprod(const Matrix& A){
	register double *newdata = new double[A.colsize * A.colsize];
	
	for (int i = 0; i < A.colsize; ++i){
	    for (int j = i; j < A.colsize; ++j){
		newdata[i*A.colsize + j] = 0.0;
		for (int k = 0; k < A.rowsize; ++k){
		    newdata[i * A.colsize + j] += A.data[k * A.colsize + i] *
			A.data[k * A.colsize + j];
		    newdata[j*A.colsize + i] = newdata[i*A.colsize + j];
		}
	    }
	}
	
	Matrix temp = Matrix (newdata, A.colsize, A.colsize);
	delete [] newdata;
	return temp;
    }
    
    // better loop ordering
    Matrix crossprod2(const Matrix& A){
	register double *newdata = new double[A.colsize * A.colsize];
	const int nr = A.rowsize;
	const int nc = A.colsize;
	
	for (int k = 0; k < nr; ++k){
	    for (int i = 0; i < nc; ++i){
		for (int j = i; j < nc; ++j){
		    newdata[i * nc + j] += A.data[k * nc + i] *
			A.data[k * nc + j];
		    newdata[j*nc + i] = newdata[i*nc + j];
		}
	    }
	}
	
	/*
	  for (int i = 0; i < A.colsize; ++i){
	  for (int j = i; j < A.colsize; ++j){
	  newdata[i*A.colsize + j] = 0.0;
	  for (int k = 0; k < A.rowsize; ++k){
	  newdata[i * A.colsize + j] += A.data[k * A.colsize + i] *
	  A.data[k * A.colsize + j];
	  newdata[j*A.colsize + i] = newdata[i*A.colsize + j];
	  }
	  }
	  }
	*/
	
	Matrix temp = Matrix (newdata, A.colsize, A.colsize);
	delete [] newdata;
	return temp;
    }
    
    
/************************    OPERATORS     ****************************/
    
// ADDITION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Addition  (Matrix + Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the sum of the two Matrices.
 */
    Matrix operator + (const Matrix & A, const Matrix & B)
    {
	if (A.rowsize == 1 && A.colsize == 1) {
	    // Case 1: A is 1 x 1 and B is n x k
	    double *sum = new double[B.size];
	    for (int i = 0; i < B.size; ++i) {
		sum[i] = A.data[0] + B.data[i];
	    }
	    Matrix temp = Matrix (sum, B.rowsize, B.colsize);
	    delete[]sum;
	    return temp;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    double *sum = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		sum[i] = A.data[i] + B.data[0];
	    }
	    Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	    delete[]sum;
	    return temp;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	    cerr << "ERROR 0055: Matrices not conformable for addition \n"
		 << "exiting program due to error" 
		 << endl;
	    exit (1);
	} else {
	    // Case 4: A is n x k and B is also n x k
	    double *sum = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		sum[i] = A.data[i] + B.data[i];
	    }
	    Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	    delete[]sum;
	    return temp;
	}
    }
    
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  Matrix + scalar
/*!
  \overload Matrix operator + (const Matrix & A, const double &b)
*/
    Matrix operator + (const Matrix & A, const double &b)
    {
	double *sum = new double[A.size];
	for (int i = 0; i < A.size; ++i) {
	    sum[i] = A.data[i] + b;
	}
	Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	delete[]sum;
	return temp;
    }
    
//  OPERATOR: Addition
//  NOTE: This operator is overloaded
//  scalar + Matrix
/*!
  \overload Matrix operator + (const double &a, const Matrix & B)
*/
    Matrix operator + (const double &a, const Matrix & B)
    {
	double *sum = new double[B.size];
	for (int i = 0; i < B.size; ++i) {
	    sum[i] = a + B.data[i];
	}
	Matrix temp = Matrix (sum, B.rowsize, B.colsize);
	delete[]sum;
	return temp;
    }
    
// SUBTRACTION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Subtraction  (Matrix - Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the difference of the two Matrices.
 */
    Matrix operator - (const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	    double *sum = new double[B.size];
	    for (int i = 0; i < B.size; ++i) {
		sum[i] = A.data[0] - B.data[i];
	    }
	    Matrix temp = Matrix (sum, B.rowsize, B.colsize);
	    delete[]sum;
	    return temp;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    double *sum = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		sum[i] = A.data[i] - B.data[0];
	    }
	    Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	    delete[]sum;
	    return temp;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	    cerr << "ERROR 0056: Matrices not conformable for subtraction\n"
		 << "exiting program due to error" 
		 << endl;
	    exit (56);
	} else {
	    // Case 4: A is n x k and B is also n x k
	    double *sum = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		sum[i] = A.data[i] - B.data[i];
	    }
	    Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	    delete[]sum;
	    return temp;;
	}
    }
    
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  Matrix - scalar
/*!
  \overload Matrix operator - (const Matrix & A, const double &b)
*/
    Matrix operator - (const Matrix & A, const double &b)
    {
	double *sum = new double[A.size];
	for (int i = 0; i < A.size; ++i) {
	    sum[i] = A.data[i] - b;
	}
	Matrix temp = Matrix (sum, A.rowsize, A.colsize);
	delete[]sum;
	return temp;
    }
    
//  OPERATOR: Subtraction
//  NOTE: This operator is overloaded
//  scalar - Matrix
/*!
  \overload Matrix operator - (const double &a, const Matrix & B)
*/
    Matrix operator - (const double &a, const Matrix & B)
    {
	double *sum = new double[B.size];
	for (int i = 0; i < B.size; ++i) {
	    sum[i] = a - B.data[i];
	}
	Matrix temp = Matrix (sum, B.rowsize, B.colsize);
	delete[]sum;
	return temp;
    }
    
// MULTIPLICATION OPERATORS IN ALL THEIR MANY FLAVORS
    
/*!
 * \brief OPERATOR: Multiplication  (Matrix * Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the product of the two Matrices.
 */
    Matrix operator *(const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	    double *prod = new double[B.size];
	    for (int i = 0; i < B.size; ++i) {
		prod[i] = A.data[0] * B.data[i];
	    }
	    Matrix temp = Matrix (prod, B.rowsize, B.colsize);
	    delete[]prod;
	    return temp;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    double *prod = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		prod[i] = A.data[i] * B.data[0];
	    }
	    Matrix temp = Matrix (prod, A.rowsize, A.colsize);
	    delete[]prod;
	    return temp;
	} else if (A.colsize != B.rowsize) {
	    // Case 3: A is n x k and B is m x j (m !=j)
	    cerr << "ERROR 0057: Matrices not conformable for multiplication\n"
		 << "exiting program due to error" 
		 << endl;
	    exit (57);
	} else {
	    // Case 4: A is n x k and B is k x j
	    register double *newdata = new double[A.rowsize * B.colsize];
	    for (int i = 0; i < A.rowsize; ++i) {
		for (int j = 0; j < B.colsize; ++j) {
		    newdata[i*B.colsize + j] = 0.0;
		    for (int k = 0; k < B.rowsize; ++k) {
			newdata[i * B.colsize + j] += A.data[i * A.colsize + k] *
			    B.data[k * B.colsize + j];
		    }
		}
	    }
	    
	    Matrix temp = Matrix (newdata, A.rowsize, B.colsize);
	    delete[]newdata;
	    return temp;
	}
    }
    
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
// Matrix * scalar
/*!
  \overload Matrix operator *(const Matrix & A, const double &b)
*/
    Matrix operator *(const Matrix & A, const double &b)
    {
	double *prod = new double[A.size];
	for (int i = 0; i < A.size; ++i) {
	    prod[i] = A.data[i] * b;
	}
	Matrix temp = Matrix (prod, A.rowsize, A.colsize);
	delete[]prod;
	return temp;
    }
    
//  OPERATOR: Multiplication
//  NOTE: This operator is overloaded
// scalar * Matrix
/*!
  \overload Matrix operator *(const double &a, const Matrix & B)
*/
    Matrix operator *(const double &a, const Matrix & B)
    {
	double *prod = new double[B.size];
	for (int i = 0; i < B.size; ++i) {
	    prod[i] = a * B.data[i];
	}
	Matrix temp = Matrix (prod, B.rowsize, B.colsize);
	delete[]prod;
	return temp;
    }
    
// ELEMENT BY ELEMENT DIVISION
/*!
 * \brief OPERATOR: Division  (Matrix / Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the dividend of each element of the two Matrices.
 */
    Matrix operator / (const Matrix & A, const Matrix & B)
    {
	// Case 1: A is 1 x 1 and B is n x k
	if (A.rowsize == 1 && A.colsize == 1) {
	    double *quot = new double[B.size];
	    for (int i = 0; i < B.size; ++i) {
		quot[i] = A.data[0] / B.data[i];
	    }
	    Matrix temp = Matrix (quot, B.rowsize, B.colsize);
	    delete[]quot;
	    return temp;
	} else if (B.rowsize == 1 && B.colsize == 1) {
	    // Case 2: A is n x k and B is 1 x 1
	    double *quot = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		quot[i] = A.data[i] / B.data[0];
	    }
	    Matrix temp = Matrix (quot, A.rowsize, A.colsize);
	    delete[]quot;
	    return temp;
	} else if (A.rowsize != B.rowsize || A.colsize != B.colsize) {
	    // Case 3: A is n x k and B is m x j (n != m or k != m)
	    cerr << "ERROR 0058: Matrices not conformable for division\n" 
		 << "exiting program due to error" 
		 << endl;
	    exit (58);
	} else {
	    // Case 4: A is n x k and B is also n x k
	    double *quot = new double[A.size];
	    for (int i = 0; i < A.size; ++i) {
		quot[i] = A.data[i] / B.data[i];
	    }
	    Matrix temp = Matrix (quot, A.rowsize, A.colsize);
	    delete[]quot;
	    return temp;;
	}
    }
    
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  Matrix / scalar
/*!
  \overload Matrix operator / (const Matrix & A, const double &b)
*/
    Matrix operator / (const Matrix & A, const double &b)
    {
	
	double *quot = new double[A.size];
	for (int i = 0; i < A.size; ++i) {
	    quot[i] = A.data[i] / b;
	}
	Matrix temp = Matrix (quot, A.rowsize, A.colsize);
	delete[]quot;
	return temp;
    }
    
//  OPERATOR: Division
//  NOTE: This operator is overloaded
//  scalar / Matrix
/*!
  \overload Matrix operator / (const double &a, const Matrix & B)
*/
    Matrix operator / (const double &a, const Matrix & B)
    {
	
	double *quot = new double[B.size];
	for (int i = 0; i < B.size; ++i) {
	    quot[i] = a / B.data[i];
	}
	Matrix temp = Matrix (quot, B.rowsize, B.colsize);
	delete[]quot;
	return temp;
    }
    
    
//!  OPERATOR: Kronecker multiplication
/*!
 * OPERATOR: Kronecker Multiplication  (Matrix % Matrix).
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return the result of the Kronecker Multiplication of the 2 Matrices.
 */
    Matrix operator % (const Matrix& A, const Matrix& B){
	double* newdata = new double[A.size * B.size];
	int count = 0;
	for (int i=0; i<A.rowsize; ++i){
	    for (int j=0; j<B.rowsize; ++j){
		for (int k=0; k<A.colsize; ++k){
		    for (int m=0; m<B.colsize; ++m){
			newdata[count] = A.data[i*A.colsize + k] * B.data[j*B.colsize + m];
			++count;
		    }
		}
	    }
	}
	Matrix temp = Matrix (newdata, A.rowsize*B.rowsize, A.colsize*B.colsize);
	delete[] newdata;
	return temp;
    }
    
    
//!  OPERATOR: Equality
/*!
 * OPERATOR: Equality.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return an integer, 1 if the Matrices are the same, 0 otherwise.
 */
    int operator == (const Matrix& A, const Matrix & B){
	if (A.rowsize != B.rowsize || A.colsize != B.colsize) return 0;
	for(int i=0; i<A.size; ++i){
	    if (A.data[i] != B.data[i]) return 0;
	}
	return 1;
    }
    
//  OPERATOR: Inequality
/*!
 * OPERATOR: Inequality
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return an integer, 1 if the Matrices are different, 0 if they are identical.
 */
    int operator != (const Matrix& A, const Matrix & B){
	return !(A==B);
    }
    
//!  OPERATOR: Element-by-element Greater Than
/*!
 * OPERATOR: Element-by-element Greater Than.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the element is 
 * greater than the other element, 0 otherwise.
 */
    Matrix operator >> (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	    cerr << "ERROR 0059: Matrices not conformable for >> operator\n" 
		 << "exiting program due to error" 
		 << endl;
	    exit(59);
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
	    double *newdata = new double[A.size];
	    for (int i = 0; i<A.size; ++i){
		newdata[i] = A.data[i] > B.data[i];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[]  newdata;
	    return temp;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] > B.data[i];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] > B.data[j];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (B.size == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.size; ++i){
		newdata[i] = A.data[i] > B.data[0];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;  
	} else {
	    cerr << "ERROR 0060: Matrices not conformable for >> operator\n"
		 << "exiting program due to error" 
		 << endl;
	    exit(60);
	}  
    }
    
//  OPERATOR: Element-by-element Greater Than
/*!
  \overload Matrix operator >> (const Matrix& A, const double& b)
*/
    Matrix operator >> (const Matrix& A, const double& b){
	double *newdata = new double[A.size];
	for (int i=0; i<A.size; ++i){
	    newdata[i] = A.data[i] > b;
	}
	Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	delete[] newdata;
	return temp;
    }
    
//!  OPERATOR: Element-by-scalar Less Than
/*!
 * OPERATOR: Element-by-element Less Than.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the element 
 * is less than the other element, 0 otherwise.
 */
    Matrix operator << (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	    cerr << "ERROR 0061: Matrices not conformable for << operator\n" 
		 << "exiting program due to error" 
		 << endl;
	    exit(61);
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
	    double *newdata = new double[A.size];
	    for (int i = 0; i<A.size; ++i){
		newdata[i] = A.data[i] < B.data[i];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[]  newdata;
	    return temp;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] < B.data[i];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] < B.data[j];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (B.size == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.size; ++i){
		newdata[i] = A.data[i] < B.data[0];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;  
	}
	
	else {
	    cerr << "ERROR 0062: Matrices not conformable for << operator" 
		 << "exiting program due to error" 
		 << endl;
	    exit(62);
	}  
    }
    
//  OPERATOR: Element-by-scalar Less Than
/*!
  \overload Matrix operator << (const Matrix& A, const double& b)
*/
    Matrix operator << (const Matrix& A, const double& b)
    {
	double *newdata = new double[A.size];
	for (int i=0; i<A.size; ++i){
	    newdata[i] = A.data[i] < b;
	}
	Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	delete[] newdata;
	return temp;
    }
    
    
//!  OPERATOR: Element-by-element Equality
/*!
 * OPERATOR: Element-by-element Equality.
 * \param A a constant reference to a Matrix \a A.
 * \param B a constant reference to a Matrix \a B.
 * \return A Matrix of 1's and 0's, 1's if the elements are equal, 
 * 0 otherwise.
 */
    Matrix operator ^= (const Matrix& A, const Matrix& B){
	if (A.rowsize != B.rowsize && A.colsize != B.colsize && (B.size > 1)){
	    cerr << "ERROR 0063: Matrices not conformable for ^= operator" 
		 << "exiting program due to error" 
		 << endl;
	    exit(63);
	}
	
	if (A.rowsize == B.rowsize && A.colsize == B.colsize){
	    double *newdata = new double[A.size];
	    for (int i = 0; i<A.size; ++i){
		newdata[i] = A.data[i] == B.data[i];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[]  newdata;
	    return temp;
	}
	
	if (A.rowsize == B.rowsize && B.colsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] == B.data[i];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (A.colsize == B.colsize && B.rowsize == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.rowsize; ++i){
		for (int j=0; j<A.colsize; ++j){
		    newdata[i*A.colsize + j] = A.data[i*A.colsize + j] == B.data[j];
		}
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;
	}
	
	if (B.size == 1){
	    double *newdata = new double[A.size];
	    for (int i=0; i<A.size; ++i){
		newdata[i] = A.data[i] == B.data[0];
	    }
	    Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	    delete[] newdata;
	    return temp;  
	}
	
	else {
	    cerr << "ERROR 0064: Matrices not conformable for ^= operator" 
		 << "exiting program due to error" 
		 << endl;
	    exit(64);
	}  
    }
    
    
//  OPERATOR: Element-by-element Equality
/*!
  \overload Matrix operator ^= (const Matrix& A, const double& b)
*/
    Matrix operator ^= (const Matrix& A, const double& b)
    {
	double *newdata = new double[A.size];
	for (int i=0; i<A.size; ++i){
	    newdata[i] = A.data[i] == b;
	}
	Matrix temp = Matrix(newdata, A.rowsize, A.colsize);
	delete[] newdata;
	return temp;
    }
    
    
    
    
} // namespace dec

#endif /* SCYTHE_DOUBLE_MATRIX_CC */
