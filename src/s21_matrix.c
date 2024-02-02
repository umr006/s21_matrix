#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error_code = OK;
  if ((rows > 0) && (columns > 0)) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    for (int i = 0; i < result->rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
    }
  } else {
    error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

void s21_remove_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    free(A->matrix[i]);
  }
  free(A->matrix);
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error_code = SUCCESS;
  if (A->rows != B->rows || A->columns != B->columns || A->rows <= 0 ||
      A->columns <= 0 || B->rows <= 0 || B->columns <= 0) {
    error_code = FAILURE;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07) {
          error_code = FAILURE;
          break;
        }
      }
      if (error_code == FAILURE) break;
    }
  }
  return error_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;
  if (!A || !B || !A->matrix || !B->matrix || A->rows <= 0 || A->columns <= 0 ||
      B->rows <= 0 || B->columns <= 0) {
    error_code = INCORRECT_MATRIX;
  } else if (A->columns != B->columns || A->rows != B->rows) {
    error_code = CALCULATION_ERROR;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return error_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;
  if (!A || !B || !A->matrix || !B->matrix || A->rows <= 0 || A->columns <= 0 ||
      B->rows <= 0 || B->columns <= 0) {
    error_code = INCORRECT_MATRIX;
  } else if (A->columns != B->columns || A->rows != B->rows) {
    error_code = CALCULATION_ERROR;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return error_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error_code = OK;
  if (!A || !A->matrix || A->rows <= 0 || A->columns <= 0) {
    error_code = INCORRECT_MATRIX;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return error_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = OK;
  if (!A || !B || !A->matrix || !B->matrix || A->rows <= 0 || A->columns <= 0 ||
      B->rows <= 0 || B->columns <= 0) {
    error_code = INCORRECT_MATRIX;
  } else if (A->columns != B->rows || A->rows != B->columns) {
    error_code = CALCULATION_ERROR;
  } else {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return error_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error_code = OK;
  if (!A->matrix || A->rows <= 0 || A->columns <= 0) {
    error_code = INCORRECT_MATRIX;
  } else {
    s21_create_matrix(A->columns, A->rows, result);

    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return error_code;
}

void find_minor(int i, int j, matrix_t A, matrix_t *result) {
  s21_create_matrix(A.rows - 1, A.columns - 1, result);
  int rows = 0;
  int columns = 0;
  for (int k = 0; k < A.rows; k++) {
    if (i == k) continue;
    for (int l = 0; l < A.columns; l++) {
      if (l == j) continue;
      result->matrix[rows][columns] = A.matrix[k][l];
      columns++;
    }
    rows++;
    columns = 0;
  }
}

int s21_determinant(matrix_t *A, double *result) {
  int error_code = OK;
  if (!A || !A->matrix || !result || A->columns <= 0 || A->rows <= 0) {
    error_code = INCORRECT_MATRIX;
  } else if (A->columns != A->rows) {
    error_code = CALCULATION_ERROR;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      *result = 0;
      for (int i = 0; i < A->columns; i++) {
        matrix_t minor_matrix;
        find_minor(0, i, *A, &minor_matrix);
        double tmp_det = 0;
        s21_determinant(&minor_matrix, &tmp_det);
        s21_remove_matrix(&minor_matrix);
        double tmp_res = pow(-1, 2 + i) * A->matrix[0][i] * tmp_det;
        *result += tmp_res;
      }
    }
  }
  return error_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error_code = OK;
  if (!A || !A->matrix || !result || A->columns <= 1 || A->rows <= 1) {
    error_code = CALCULATION_ERROR;
  } else if (A->columns != A->rows) {
    error_code = CALCULATION_ERROR;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        matrix_t minor;
        find_minor(i, j, *A, &minor);
        double det = 0;
        s21_determinant(&minor, &det);
        det *= pow(-1, 2 + i + j);
        result->matrix[i][j] = det;
        s21_remove_matrix(&minor);
      }
    }
  }
  return error_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error_code = OK;
  if (!A || !A->matrix || !result || A->columns <= 1 || A->rows <= 1) {
    error_code = INCORRECT_MATRIX;
  } else if (A->rows != A->columns) {
    error_code = CALCULATION_ERROR;
  } else {
    double determ = 0;
    matrix_t tmp;
    matrix_t tmp1;
    s21_determinant(A, &determ);
    if (fabs(determ) > 1e-7) {
      s21_calc_complements(A, &tmp);
      s21_transpose(&tmp, &tmp1);
      s21_remove_matrix(&tmp);
      s21_mult_number(&tmp1, 1. / fabs(determ), result);
      s21_remove_matrix(&tmp1);
      for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->columns; j++) {
          result->matrix[i][j] *= -1;
        }
      }
    } else {
      error_code = CALCULATION_ERROR;
    }
  }
  return error_code;
}