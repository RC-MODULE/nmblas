 /**
 * \brief   nmblas (NeuroMatrix Basic Linear Algebra Subroutines)
 * \details  This class is used to demonstrate a number of section commands.
 * \author  Alexandr Bolotnikov 
 * \version  1.0
 * \date   2017-05-23T15:13:13
 * \bug    
 * \warning  
 * \copyright (c) RC Module Inc.
 * \file nmblas.h
 */
 
 
 //! \defgroup LEVEL1 BLASS-LEVEL1
 //! \{
 
 //! \}


#ifndef _NMBLAS_H_INCLUDED_
#define _NMBLAS_H_INCLUDED_

#ifdef __cplusplus
		extern "C" {
#endif	



enum nm_trans{nm_n=0,nm_t=1};

 /**
 * \ingroup LEVEL1
 *
 * \brief Функция почленно складывает элементы вектора по модулю
 * 
 * \param [in] N количество элементов в векторе 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \return сумма элементов вектора по модулю
 */
double nmblas_dasum(
	const int N,
	const double *X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief выполняет операцию Y = ALPHA*X+Y, где ALPHA скаляр, а X и Y вектора
 * 
 * \param [in] N количество элементов в векторе 
 * \param [in] ALPHA скаляр
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив 
 * \param [in] INCY смещение между элементами в массиве 
 */
void nmblas_daxpy(
	const int N,
	const double ALPHA,
	const double *X,
	const int INCX,
	double *Y,
	const int INCY
); 
/**
 * \ingroup LEVEL1
 * \brief копирует элементы массива X в массив Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY смещение между элементами в массиве 
 */

void nmblas_dcopy(
	const int N,
	const double* X,
	const int INCX,
	double* Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 * 
 * \brief возвращает скалярное произведение векторов X и Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY смещение между элементами в массиве 
 * \return скалярное произведение X и Y 
 */
double nmblas_ddot(
	const int N,
	const double* X,
	const int INCX,
	const double* Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief возвращает Евклидову норму вектора X
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 *
 * \return Евклидова норма вектора 
 */
double nmblas_dnrm2(
	const int N,
	const double *X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief выполняет поворот массива точек в двумерном пространстве (поворот Гивена)<br>
 * \brief [ x_i ] = [ c  s ] [ x_i ]<br>
 * \brief [ y_i ] = [ -s c ] [ y_i ]<br>   
 * \brief Если n <= 0 или если c = 1 и s = 0, функция выходит сразу.	
 * 
 * \param [in] N количество точек в массивах 
 * \param [in,out] *X указатель на массив, который содержит x координату точки 
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив, который содержит y координату точки 
 * \param [in] INCY смещение между элементами в массиве 
 * \param [in] С значение косинуса угла поворота
 * \param [in] S значение синуса угла поворота
 *
 */
void nmblas_drot(
	const int N,
	double *X,
	const int INCX,
	double *Y,
	const int INCY,
	const double C,
	const double S
); 
/**
 * \ingroup LEVEL1
 *
 * \brief вычисяет параметры необходимые для пороворота Гивена 
 * 
 * \param [in,out] *A указатель на параметр a, после выполнения a будет переписан параметром r(смотри детальное описание) 
 * \param [in,out] *B указатель на параметр b, после выполнения a будет переписан параметром z(смотри детальное описание)
 * \param [out] *C после выполнения будет записан параметр с(смотри делальное описание) 
 * \param [out] *S после выполнения будет записан параметр s(смотри делальное описание)
 *
 * \details  
 * drotg вычисляет параметрв для поворота Гивена.<br>
 * Принимая скаляр а и b,программа вычисляет следущие параметры:<br>
 *
 *    sigma = sgn(a) если |a| > |b|, и sgn(b) в противном случаии;<br>
 *    r  		= sigma * sqrt( a^2 + b^2 );<br>
 *    c  		= a / r если r <> 0, иначе 1;<br>
 *    s  		= b / r если r <> 0, иначе 0;<br><br>
 * \details 
 * числа c, s и r будут удовлетворять следующему матричному уровнению:<br>
 *
 *    [ c   s ] [ a ] =  [ r ]<br>
 *    [ -s c ] [ b ] =  [ 0 ]<br>
 *
 *	Sigma не важна для вычислений матрицы поворота Гивена, но важна для<br> 
 *	дальнейшего надежного востановления c и s из одного сохраненого числа.<br>
 *  Для этого данная функция также вычисляет z как показано ниже  <br>   
 *
 *    z =  s        если |a| > |b|,<br>
 *    z =  1 / c    если |b| >= |a| и c <> 0,1<br>
 *    z =  1        если c = 0.<br><br>
 * \details 
 * Функция возвращает r, переписывая a и z переписывая b, также возвращая c и s.<br>
 * Если позже будет необходимость в реконструировании c и s из z это можно будет сделать следующим образом<br>
 *
 *    если  z  = 1, установить c = 0             и s = 1,<br>
 *    если |z| < 1, установить c = sqrt(1 - z^2) и s = z,<br>
 *    если |z| > 1, установить c = 1 / z         и s = sqrt(1 - c^2).<br><br>
 * \details 
 * Смотри ``Basic Linear Algebra Subprograms for Fortran Usage''  C. Law-<br>
 * son, R. Hanson, D. Kincaid and F. Krogh, ACM Transactions on Mathema-<br>
 * tical Software, 1979, 5(3) pp 308-323, для болеее подробной информации.<br>
 */ 
void nmblas_drotg(
	double *a, 
	double *b, 
	double *c, 
	double *s
);
/**
 * \ingroup LEVEL1
 * \brief функция применяет модифицированный поворот Гивена (modified-Givens rotation), матрица H предается через *PARAM<br>
 * \brief [ x_i ]   [ h_11 h_12 ] [ x_i ]<br>
 * \brief [ y_i ] = [ h_21 h_22 ] [ y_i ]<br>   
 * \brief Если n <= 0 или если H индентичная матрица, выход из функции происходит сразу
 *
 * \param [in] N пеередает длину векторов.
 * \param [in,out] *X указатель на вектор X, после выполнения будет переписан результатом
 * \param [in] INCX передает инкремент массива X (шаг между элементами)
 * \param [in,out] *Y указатель на вектор Y, после выполнения будет переписан результатом
 * \param [in] INCY передает инкремент массива Y (шаг между элементами)
 * \param [in] *PARAM указатель на массив размерности не менее 5(double), который содержит матрицу H. 
 *
*/
void nmblas_drotm(
	const int N,
	double *X,
	const int INCX,
	double *Y,
	const int INCY,
	double *param
);
/**
 * \ingroup LEVEL1
 * 
 * \brief выполняет умножение вектора на скаляр X = alpha*X
 * 
 * \param [in] N количество точек в массивах 
 * \param [in] ALPHA значение скаляра 
 * \param [in,out] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 *
 */  
void nmblas_dscal(
	const int N,
	const double ALPHA,
	double*X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 * \brief возвращает скалярное произведение векторов X и Y тапа float. Вычисления внутри функции происходят в числах двойной(double) точности. 
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY смещение между элементами в массиве 
 * \return double скалярное произведение X и Y 
 */
double nmblas_dsdot(
	const int N,
	const float* X,
	const int INCX,
	const float* Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief выполняет перестановку массивов X и Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in,out] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY смещение между элементами в массиве 
 *	
 */
void nmblas_dswap(
	const int N,
	double *X,
	const int INCX,
	double *Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief возвращает индекс наибольшего по модулю элемента в массиве
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX смещение между элементами в массиве 
 * \return индекс элемента 
 */ 
int nmblas_idamax(
	const int N,
	const double* X,
	const int INCX
);
//////////////////////////////////////////////////////level 1 single precition
/**
 * \ingroup LEVEL1
 *
 * \brief возвращает индекс наибольшего по модулю элемента в массиве
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив 
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \return индекс элемента 
 */ 
int nmblas_isamax(
	const int N,
	const float* X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief Функция почленно складывает элементы вектора по модулю
 * 
 * \param [in] N количество элементов в векторе 
 * \param [in] *X указатель на массив 
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \return сумма элементов вектора по модулю
 */
float nmblas_sasum(
	const int N,
	const float *X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief выполняет операцию Y = ALPHA*X+Y, где ALPHA скаляр, а X и Y вектора
 * 
 * \param [in] N количество элементов в векторе 
 * \param [in] ALPHA скаляр
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY ДОЛЖЕН быть равен 1 
 */
void nmblas_saxpy(
	const int N,
	const float ALPHA,
	const float *X,
	const int INCX,
	float *Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief возвращает Евклидову норму вектора X
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 *   
 */  
float nmblas_scnrm2(
	const int N,
	const float* X,
	const int INCX
);

/**
 * \ingroup LEVEL1
 *
 * \brief копирует элементы массива X в массив Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY ДОЛЖЕН быть равен 1 
 *	
 */
void nmblas_scopy(
	const int N,
	const float* X,
	const int INCX,
	float* Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief возвращает скалярное произведение векторов X и Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \param [in] *Y указатель на массив
 * \param [in] INCY ДОЛЖЕН быть равен 1 
 * \return скалярное произведение X и Y 
 */
float nmblas_sdot(
	const int n,
	const float *x,
	const int inc_x,
	const float *y,
	const int inc_y
);
/**
 * \ingroup LEVEL1
 *
 * \brief Bозвращает скалярное произведение векторов X и Y с добавлением константы B.<br>
 * \brief Внутри функции используется двойная(double) точность
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY ДОЛЖЕН быть равен 1 
 * \return скалярное произведение X и Y + B 
 */
float nmblas_sdsdot(
	const int N,
	const float B,
	const float* X,
	const int INCX,
	const float* Y,
	const int INCY
);
/**
 * \ingroup LEVEL1
 *
 * \brief Bозвращает Евклидову норму вектора X
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равен 1
 * \return Bозвращает Евклидову норму вектора X 
*/
float nmblas_snrm2(
	const int N,
	const float* X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief Bыполняет поворот массива точек в двумерном пространстве
 * 
 * \param [in] N количество точек в массивах 
 * \param [in,out] *X указатель на массив, который содержит x координату точки 
 * \param [in] INCX ДОЛЖЕН быть равен 1 
 * \param [in,out] *Y указатель на массив, который содержит y координату точки 
 * \param [in] INCY ДОЛЖЕН быть равен 1 
 * \param [in] С значение косинуса угла поворота
 * \param [in] S значение синуса угла поворота
 */
void nmblas_srot(
	const int N,
	float *X,
	const int INCX,
	float *Y,
	const int INCY,
	const float C,
	const float S
);
/**
 * \ingroup LEVEL1
 *
 * \brief функция применяет модифицированный поворот Гивена (modified-Givens rotation),<br> 
 * \brief матрица H предается через *PARAM<br><br>
 * \brief  [ x_i ] = [ h_11 h_12 ] [ x_i ]<br>
 * \brief  [ y_i ] = [ h_21 h_22 ] [ y_i ]<br> 
 * \brief Если n <= 0 или если H индентичная матрица, выход из функции происходит сразу
 *
 * \param [in] N пеередает длину векторов.
 * \param [in,out] *X указатель на вектор X, после выполнения будет переписан результатом
 * \param [in] INCX ДОЛЖЕН быть равным 1
 * \param [in,out] *Y указатель на вектор Y, после выполнения будет переписан результатом
 * \param [in] INCY ДОЛЖЕН быть равным 1
 * \param [in] *PARAM указатель на массив размерности не менее 5(float), который содержит матрицу H. 
*/
void nmblas_srotm(
	const int N,
	float *X,
	const int INCX,
	float *Y,
	const int INCY,
	float *param
);
 /**
 * \ingroup LEVEL1
 *
 * \brief выполняет умножение вектора на скаляр X = ALPHA*X
 * 
 * \param [in] N количество точек в массивах 
 * \param [in] ALPHA значение скаляра 
 * \param [in] *X указатель на массив, который содержит x координату точки 
 * \param [in] INCX ДОЛЖЕН быть равным 1 
 */  
void nmblas_sscal(
	const int N,
	const float ALPHA,
	float *X,
	const int INCX
);
/**
 * \ingroup LEVEL1
 *
 * \brief выполняет перестановку массивов X и Y
 * 
 * \param [in] N количество элементов в массиве 
 * \param [in] *X указатель на массив
 * \param [in] INCX ДОЛЖЕН быть равным 1
 * \param [in,out] *Y указатель на массив
 * \param [in] INCY ДОЛЖЕН быть равным 1 
 */
void nmblas_sswap(
	const int N,
	const float *X,
	const int INCX,
	const float *Y,
	const int INCY
);

/////////////////LEVEL 2 
/**
 * \ingroup LEVEL 2
 *
 * \brief выполняет следующие операцию Y = ALPHA*A[M,N]*X + BETTA*Y,<br>
 * \brief где A матрица размера M на N,ALPHA и BETTA скаляры, X и Y вектора корректной длины 
 * 
 * \param [in] TRANS при значение 0 матрица A загружается не транспонировано, при значение 1 матрица загружается транспонировано 
 * \param [in] M количество рядов матрицы A
 * \param [in] N количество строк матрицы A
 * \param [in] ALPHA скаляр ALPHA
 * \param [in] *A указатель на массив матрицы A 
 * \param [in] LDA stride массива матрицы 
 * \param [in] *X указатель на массив вектора X
 * \param [in] INCX ДОЛЖЕН быть равным 1 
 * \param [in] BETTA скаляр BETTA 
 * \param [in] *Y указатель на массив вектора Y
 * \param [in] INCY ДОЛЖЕН быть равным 1 
 */
void nmblas_sgemv(
  const enum nm_trans    TRANS,
  const int         M,
  const int         N,
  const float        ALPHA,
  const float        * A,
  const int         LDA,
  const float        * X,
  const int         INCX,
  const float        BETA,
	    float        * Y,
  const int         INCY
);
/**
 * \ingroup LEVEL 2
 *
 * \brief выполняет следующие операцию Y = ALPHA*A[M,N]*X + BETTA*Y,<br>
 * \brief где A матрица размера M на N,ALPHA и BETTA скаляры, X и Y вектора подходящей длины 
 * 
 * \param [in] TRANS при значение 0 матрица A загружается не транспонировано, при значение 1 матрица загружается транспонировано 
 * \param [in] M количество рядов матрицы A
 * \param [in] N количество строк матрицы A
 * \param [in] ALPHA скаляр ALPHA
 * \param [in] *A указатель на массив матрицы A 
 * \param [in] LDA stride массива матрицы 
 * \param [in] *X указатель на массив вектора X
 * \param [in] INCX ДОЛЖЕН быть равным 1 
 * \param [in] BETTA скаляр BETTA 
 * \param [in] *Y указатель на массив вектора Y
 * \param [in] INCY ДОЛЖЕН быть равным 1 
 */
void nmblas_dgemv(
  const enum nm_trans    TRANS,
  const int         M,
  const int         N,
  const double        ALPHA,
  const double        * A,
  const int         LDA,
  const double        * X,
  const int         INCX,
  const double        BETA,
        double     * Y,
  const int         INCY
);
 
////LEVEL 3
#ifdef __cplusplus
		};
#endif

void MullMatrix_f( 
	const float* 	A, 
	int 	pI, 
	int 	ldA, 
	const float* 	B,
  int 	pK, 
  int 	ldB,
  float *C,
  int 	pJ, 
  int 	ldC,
  bool 	plusC 
);


	
#endif // _INIT_H_INCLUDED_