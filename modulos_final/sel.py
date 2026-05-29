import numpy as np
import time

def eliminacion_DD(A, b):
    matrix_a = np.insert(A, A.shape[0], b, 1)
    x_solucion = np.zeros_like(b)

    n = len(b)

    for j in range(n + 1):
        for i in range(j + 1, n):
            factor = matrix_a[i, j] / matrix_a[j, j]
            matrix_a[i, 0:n+1] = matrix_a[i, 0:n+1] - factor * matrix_a[j, 0:n+1]

    for k in range(n-1, -1, -1):
        x_solucion[k] = (matrix_a[k, n] - np.dot(matrix_a[k, k+1:n], x_solucion[k+1:n])) / matrix_a[k, k]

    return x_solucion


def Eliminacion_con_pivote(A,b):
    matriz_a = np.insert(A, A.shape[0], b, 1)
    x_solucion = np.zeros_like(b)
    n = len(b)
    for j in range(n):
        if matriz_a[j,j] == 0:
            copia_actual = matriz_a[j:n+1]
            print(matriz_a[j])
            matriz_a[j] = matriz_a[0:n+1]
            matriz_a[j+1]= copia_actual

        for i in range(j+1,n):
            factor = matriz_a[i,j]/matriz_a[j,j]
            matriz_a[i, 0:n+1] = matriz_a[i, 0:n+1] - factor*matriz_a[j, 0:n+1]

    for k in range(n-1,- 1,-1):
        x_solucion[k] = (matriz_a[k,n]-np.dot(matriz_a[k,k+1:n],x_solucion[k+1:n]))/matriz_a[k,k]
    return x_solucion

def jacobi_matrices(A,b, tol):
  x0=np.array([1,0,1,0],float)
  D=np.diag(np.diag(A))
  U= D-np.triu(A)
  L=D-np.tril(A)

  T_jab=np.dot(np.linalg.inv(D),L+U)
  c_jab=np.dot(np.linalg.inv(D),b)

  eigvalues, eigvectores = np.linalg.eig(T_jab)
  radio_espectral=max(abs(eigvalues))

  if (radio_espectral >= 1):
    return
  else:
    error = 1
    cont=0
    tiempo_inicial=time.time()
    while error > tol:
      cont+=1
      x1 = np.dot(T_jab,x0)+c_jab
      error = max(abs(x1-x0))
      if cont == 3:
        x0=x1
    tiempo_final=time.time()
    tiempo_total=tiempo_final-tiempo_inicial
    return x1, cont, tiempo_total



def jacobi_sumas(A,b,x0,tol):
  Nmax=50
  conteo=0
  error = 1
  x_new=np.zeros_like(b)
  while error > tol and conteo<Nmax:
    for i in range(len(b)):
      suma=0
      for j in range(len(b)):
        if j!=i:
          suma+=A[i,j]*x0[j]
      x_new[i] = (b[i]-suma) / A[i,i]
    conteo+=1
    error = max(abs(x_new-x0))
    x0=x_new.copy()
  return x_new, conteo


def gauss_seidel_matrices(A, b, x0, tol):
    D = np.diag(np.diag(A))
    U = D - np.triu(A)
    L = D - np.tril(A)
    T_g = np.dot(np.linalg.inv(D - L), U)
    C_g = np.dot(np.linalg.inv(D - L), b)

    eigvalues, eigvectors = np.linalg.eig(T_g)
    radio_espectral = np.max(np.abs(eigvalues))

    if radio_espectral >= 1:
        return None, None
    error = 1
    iteraciones = 0
    while error > tol:
        iteraciones += 1
        x_new = np.dot(T_g, x0) + C_g
        error = np.max(np.abs(x_new - x0))
        x0 = x_new
    return x_new, iteraciones


def Gauss_seidel_sumas(A,b,x0,tol):
  n=len(b)
  error = 1
  while error > tol:
    x_new=np.zeros_like(x0)
    for i in range(n):
      sum1=0
      for j in range(i):
        sum1+=A[i,j]*x_new[j]
      sum2=0
      for j in range(i+1,n):
        sum2+=A[i,j]*x0[j]
      x_new[i]= (b[i]-sum1-sum2)/A[i,i]
    error = max(abs(x_new-x0))
    x0=x_new
  return x_new