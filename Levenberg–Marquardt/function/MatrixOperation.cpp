#include "MatrixOperation.h"
#include <stdio.h>
#include <memory>

void MatrixOperation::MatrixPrint(double** matrix, int row, int col)
{
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%f, ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double** MatrixOperation::Transpose(double** matrix, double** buffer, int row, int col)
{
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            buffer[i][j] = matrix[j][i];
        }
    }
    return buffer;
}

double** MatrixOperation::Multiple(double** A, double** B, double** buffer, int A_row, int B_row, int B_col)
{

    if (B_col != 1)
    {
        for (int i = 0; i < A_row; i++)
        {
            for (int j = 0; j < B_col; j++)
            {
                buffer[i][j] = 0;
                for (int k = 0; k < B_row; k++)
                {
                    (buffer[i][j]) += (A[i][k] * B[k][j]);
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < A_row; i++)
        {
            for (int j = 0; j < B_col; j++)
            {
                buffer[j][i] = 0;
                for (int k = 0; k < B_row; k++)
                {
                    (buffer[j][i]) += (A[i][k] * B[j][k]);
                }
            }
        }
    }

    return buffer;
}

void MatrixOperation::VectorDot(double* A, double* B, double* buffer,int length)
{
    buffer[0] = 0;
    for (int i = 0; i < length; i++)
    {
        *buffer += A[i] * B[i];
    }
}

double** MatrixOperation::UnitMartix(int row, int col, double value)
{
    double** _array = (double**)calloc(row, sizeof(double*));
    for (int i = 0; i < row; i++) {
        _array[i] = (double*)calloc(col, sizeof(double));
    }

    for (int i = 0; i < row; i++)
    {
        //for (int j = 0; j < col; j++)
        //{
            _array[i][i] = value;
        //}
    }


    return _array;
}

void MatrixOperation::MatrixAdd(double** augend,double** A, double** B,int row,int col)
{
    for (uint16_t i = 0; i < row; i++)
    {
        for (uint16_t j = 0; j < col; j++)
           A[i][j] = augend[i][j] + B[i][j];
    }


}




double** MatrixOperation::InverseMatrix(double** A,int element)
{
    double tmp;
    

    double** Output = (double**)calloc(element,sizeof(double*));
    for (int i = 0; i < element; i++) {
        Output[i] = (double*)calloc(element,sizeof(double));
    }

    for (int i = 0; i < element; i++)
    {
        Output[i][i] = 1;
    }

    /* 1 0 0 */
    /* 0 1 0 */
    /* 0 0 1 */ /* 高斯消去法 */

    int x, y, k;

    for (k = 0; k < element; k++)
    {
        if (A[k][k] == 0)
        {
            for (x = k + 1; x < element; x++)  //判斷是否對角元素是否為0
            {
                for (y = 0; y < element; y++) //swap chang line
                {
                    tmp = A[k][y];
                    A[k][y] = A[x][y];
                    A[x][y] = tmp;
                    tmp = Output[k][y];
                    Output[k][y] = Output[x][y];
                    Output[x][y] = tmp;
                }

                if (A[k][k] != 0)
                    break;
            }

            if (A[k][k] == 0)
            {
                // 無反矩陣
                return nullptr;
            }
        }


        tmp = 1 / A[k][k];
        for (y = 0; y < element; y++)
        {
            Output[k][y] *= tmp;
            A[k][y] *= tmp;
        }

        for (x = 0; x < element; x++)
        {
            if (x == k)
                continue;

            tmp = A[x][k] / A[k][k]; // A[x][k]
            for (y = 0; y < element; y++)
            {
                A[x][y] -= tmp * A[k][y];
                Output[x][y] -= tmp * Output[k][y];
            }
        }

        tmp = 1 / A[k][k];

        for (y = 0; y < element; y++)
        {
            A[k][y] *= tmp;
            Output[k][y] *= tmp;
        }
        //MatrixPrint(Output, 4, 4);
    }
     return Output;

}
