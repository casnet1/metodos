#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<conio.h>
#include<string.h>
#define PI 3.14159265358979323846

#define RAIZ(x) sin(x)*cos(x)
#define RAIZ1(x) exp(x)*cos(x)
#define RAIZ2(x) pow(x,2)*pow(cos(x),2)
#define RAIZ3(x) pow(sin(x),3)*(x-5)*(x+3)
/*
#define RAIZ(x) (pow(x,2)+2*x-15)*exp(x)
#define RAIZ1(x) pow(x,4)-8*pow(x,3)-35*pow(x,2)-450*x-1001
#define RAIZ2(x) pow(sin(3*x-PI),4)
#define RAIZ4(x)   0.5*pow(x,3)-4*pow(x,2)+5.5*x-1
*/
#define DER(x)     1.5*pow(x,2)-8*x+5.5


float** creaMatriz(int n){//Funcion que reserva la memoria para la matriz
	int i;
	float** A;
	
	do{	
		A = (float**)malloc(n*sizeof(float*));
	}while(A == NULL);
	
	for(i = 0; i < n; i++){
		do{
			A[i] = (float*)calloc(n,sizeof(float));
		}while(A[i] == NULL);
	}
	
	return A;
}

void ingresaMatriz(float** A,int n){//Funcion que hace la lectura de la matriz 
	int i,j;
	
	system("cls");
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("A[%i][%i] : ",i+1,j+1);
			scanf("%f",&A[i][j]);
		}
	}
	
	printf("\n\n");
	
	system("pause");
	system("cls");
}

float** igualar(float** A,int n){//Funcion que iguala dos matrices
	float** B;
	int i,j;
	
	do{
		B = (float**)malloc(n*sizeof(float*));
	}while(B == NULL);
	
	for(i = 0; i < n; i++){
		B[i] = (float*)malloc(n*sizeof(float));
	}
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			B[i][j] = A[i][j];
		}
	}
	
	return B;
}

int cero(float *A,int n){//Booleana que determina si hay puros ceros en un vector
	int cont = 0,i;
	
	for(i = 0; i < n; i++){
		if(A[i] == 0)
			cont++;
	}
	
	if(cont == n)
		return 1; //Todos los elementos de ese vector son 0
	else
		return 0; //No todos los elementos o ninguno de ese vector son 0
}

int equal(float* A, float *B,int n){//Booleana que determina si dos matrices son iguales
	int flag = 1,i;
	
	for(i = 0; i < n; i++){
		if(A[i] != B[i]){
			flag = 0;
			break;
		}
	}
	
	return flag;
}

void intercambio(float **B, int r1, int r2){//Funcion que intercambia renglones
	float *aux;
	
	aux = B[r1];
	B[r2] = B[r1];
	B[r1] = aux;
	
	
}

void suma(float **A,float esc,int f, int d,int n){//Funcion que hace la suma de dos renglones
	float *C;
	int i;
	
	C = (float*)malloc(n*sizeof(float));
	
	for(i = 0; i < n; i++){
		C[i] = A[f][i];
		C[i] *= esc;
	}
	for(i = 0; i < n; i++){
		A[d][i] += C[i];
	}
	
	
	free(C);
}

void lib(float **A,int n){//Funcion que libera la memoria de la matriz
	int i;
	
	for(i = 0; i < n; i++){
		free(A[i]);
	}
	
	free(A);
}

float det(float** A,int n){//Funcion que retorna el determinante de la matriz
	float **B;
	float d = 1;
	float alfa;
	int j,i;
	
	B = igualar(A,n);
	
	for(j = 0; j < n-1; j++){
		if(cero(B[j],n) || equal(B[j],B[j+1],n)){
			d = 0;
			break;
		}
		else if(B[j][j] == 0){
			d *= -1;
			intercambio(B,j,j+1);
			j--;
		}
		else{
			for(i = j; i < n; i++){
				if(i != j){
					alfa = (-1)*(B[i][j])/B[j][j];
					suma(B,alfa,j,i,n);
				}
			}
			d *= B[j][j];
		}
	}
	d*= B[n-1][n-1];
	lib(B,n);
	
	return d;
}

void imp_mat(float** A,int n){//Funcion que imprime la matriz
	int i,j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("[%.5f] ",A[i][j]);
		}
		printf("\n");
	}
	
	printf("\n\n");

}

int sub_det(float** A, int n){//Booleana que determina si la matriz es definida positiva
	int i,k,l,flag = 1;
	float** C;
	
	for(i = 0; i < n && flag; i++){
		C = creaMatriz(i+1);
		for(k = 0; k < i+1; k++){
			for(l = 0; l < i+1; l++){
				C[k][l] = A[k][l];
			}
		}
		if(det(C,i+1) <= 0)
			flag = 0;
		lib(C,i+1);
	}
	
	return flag;
}

int simetrica(float** A, int n){//Booleana que determina si la matriz es simetrica
	int i,j;
	int flag = 1;
	
	for(i = 0; i < n && flag; i++){
		for(j = i; j < n; j++){
			if(j != i){
				if(A[i][j] != A[j][i])
					flag = 0;
			}
		}
	}
		
	return flag;
}

float** identidad(int n){//Funcion que asigna a una matriz la identidad de nxn
	int i,j;
	float** I;
	
	do{
		I = (float**)malloc(n*sizeof(float*));
	}while(I == NULL);
	
	for(i = 0; i < n; i++){
		do{
			I[i] = (float*)calloc(n,sizeof(float));
		}while(I[i] == NULL);
		I[i][i] = 1;
	}
	
	return I;
}

float sum_cuad(float** L,int i){//Funcion que retorna la suma de los cuadrados necesaria para los valores de L
	float r = 0;
	int j;
	
	for(j = 0; j<= i-1; j++){
		r += pow(L[i][j],2);
	}
	
	return r;
}

float suma_i_j_esima(float** L, int i, int j){//Funcion que calcula la suma necesaria para los valores de L
	float r = 0;
	int k;

	for(k = 0; k <= i-1; k++){
		r += L[i][k]*L[j][k];
	}
	return r;	
}
 	
void diagonal(float** L,float** A,int i,int n){//Funcion que calcula el elemento de la diagonal del renglon deseado
	L[i][i] = sqrt(A[i][i]-sum_cuad(L,i));
}

void elemento_iesimo(float** L, float** A, int i, int n){//Funcion que calcula cada elemento iesimo de L
	int j;
	
	for(j = i+1; j < n; j++){
		L[j][i] = (1/L[i][i])*(A[j][i] - suma_i_j_esima(L,j,i));
	}
}

void transponer(float** LT, float** L, int n){//Funcion que asgina a una matriz la transpuesta de otra
	int i,j;
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			LT[i][j] = L[j][i];
		}
	}
	
}

int cholesky (){
	float** A;
	float** L;
	float** LT;
	int n;
	int i;
	
	printf("Introduce el orden de la Matriz: ");
	scanf("%i",&n);
	n = 5;
	
	A = creaMatriz(n);
	ingresaMatriz(A,n);
	
	if(!sub_det(A,n)){
		printf("La matriz ingresada no puede ser factorizada por Cholseky ya que no es definida positiva.\n\n");
		system("pause");
		system("cls");
		lib(A,n);
	}
	else if(!simetrica(A,n)){
		printf("La matriz si es definida positiva pero no cumple con ser simetrica para este metodo.\n\n");
		system("pause");
		system("cls");
		lib(A,n);
	}
	else{
		L = identidad(n);
		LT = identidad(n);
		
		L[0][0] = sqrt(A[0][0]);
		
		for(i = 1; i < n; i++){
			L[i][0] = A[i][0]/L[0][0];
		}
		
		for(i = 1; i < n; i++){
			diagonal(L,A,i,n);
			elemento_iesimo(L,A,i,n);
		}
		
		transponer(LT,L,n);
		
		system("cls");
		
		printf("Matriz Original: \n\n");
		imp_mat(A,n);
		
		printf("Matriz L:\n\n");
		imp_mat(L,n);
		
		printf("Matriz U: \n\n");
		imp_mat(LT,n);
		
		lib(A,n);
		lib(L,n);
		lib(LT,n);
		
		system("pause");
	}
	
}

float** triang(int n){//Funcion que vuelve cero toda una matriz
	int i,j;
	float** I;
	
	do{
		I = (float**)malloc(n*sizeof(float*));
	}while(I == NULL);
	
	for(i = 0; i < n; i++){
		do{
			I[i] = (float*)calloc(n,sizeof(float));
		}while(I[i] == NULL);
	}
	
	return I;
}

void Dool(float** L,int n){//Funcion especial del metodo de Doolittle que asigna un 1 a toda la diagonal de L
	int i;
	
	for(i = 0; i < n; i++){
		L[i][i] = 1;
	}
}

void primer_renglon(float** U, float** A, int n){//Funcion que calcula los valores del primer renglon de U
	int i;
	
	for(i = 0; i < n; i++){
		U[0][i] = A[0][i];
	}
}

void primera_columna(float** L, float** A, float** U,int n){//Funcion que calcula la primera columna de L
	int i;
	
	for(i = 1; i < n; i++){
		L[i][0] = A[i][0]/U[0][0];
	}
}

float suma_i_j_esimad(float** L, float** U, int i, int j){//Funcion que retorna la suma necesaria para cada elemento de U y de L
	float r = 0;
	int k;

	for(k = 0; k <= i-1; k++){
		r += L[i][k]*U[k][j];
	}
	return r;	
}

void fila_iesima(float** U, float** L, float** A,int i,int n){//Funcion que calcula toda la fila iesima de U
	int j;
	
	for(j = i; j < n; j++){
		U[i][j] = A[i][j] - suma_i_j_esimad(L,U,i,j);
	}
}

void col_iesima(float** L, float** U, float** A,int i,int n){//Funcion que calcula toda la columna iesima de L
	int j;
	
	for(j = i+1;j < n; j++){
		L[j][i] = (1/U[i][i])*(A[j][i] - suma_i_j_esimad(L,U,j,i));
	}
}

int dolittle(){
	
	float** A;
	float** L;
	float** U;
	int n,i;
	
	printf("Introduce el orden de la matriz: ");
	scanf("%i",&n);
	
	A = creaMatriz(n);

	ingresaMatriz(A,n);
	
	if(!sub_det(A,n)){
		printf("No todos los determinantes de las submatrices son diferentes de cero, no se puede factorizar la matriz.\n\n");
		system("pause");
		system("cls");
		lib(A,n);
	}
	else{
		L = triang(n);
		U = triang(n);
		Dool(L,n);
		
		primer_renglon(U,A,n);
		primera_columna(L,A,U,n);

		for(i = 1; i < n; i++){
			fila_iesima(U,L,A,i,n);
			col_iesima(L,U,A,i,n);
		}
		
		printf("Matriz Original: \n\n");
		imp_mat(A,n);
		
		printf("L: \n\n");
		imp_mat(L,n);
		
		printf("U: \n\n");
		imp_mat(U,n);
		
		system("pause");
		
		lib(A,n);
		lib(L,n);
		lib(U,n);
		
	}
}

void primera_columnaCr(float** L, float** A,int n){//Funcion que calcula la primera columna de L
	int i;
	
	for(i = 0; i < n; i++){
		L[i][0] = A[i][0];
	}
}

void primer_renglonCr(float** U, float** A, float** L, int n){//Funcion que calcula el primer renglon de U
	int i;
	
	for(i = 1; i < n; i++){
		U[0][i] = A[0][i]/L[0][0];
	}
}

void Crout(float** U,int n){//Funcion especial del metodo de Crout que asigna un 1 a cada elemento de la diagonal de U
	int i;
	
	for(i = 0; i < n; i++){
		U[i][i] = 1;
	}
}

void col_iesimaCr(float** L, float** U, float** A,int i,int n){//Funcion que calcula la columna iesima de L
	int j;
	
	for(j = i;j < n; j++){
		L[j][i] = A[j][i] - suma_i_j_esimad(L,U,j,i);
	}
}

void fila_iesimaCr(float** U, float** L, float** A,int i,int n){//Funcion que calcula la fila iesima de U
	int j;
	
	for(j = i+1; j < n; j++){
		U[i][j] = (1/L[i][i])*(A[i][j] - suma_i_j_esimad(L,U,i,j));
	}
}

int Crou(){
	float** A;
	float** L;
	float** U;
	int n,i;
	
	printf("Introduce el orden de la matriz: ");
	scanf("%i",&n);

	A = creaMatriz(n);

	ingresaMatriz(A,n);
	
	if(!sub_det(A,n)){
		printf("No todos los determinantes de las submatrices son diferentes de cero, no se puede factorizar la matriz.\n\n");
		system("pause");
		system("cls");
		lib(A,n);
	}
	else{
		L = triang(n);
		U = triang(n);
		Crout(U,n);
		
		primera_columnaCr(L,A,n);
		primer_renglonCr(U,A,L,n);

		for(i = 1; i < n; i++){
			col_iesimaCr(L,U,A,i,n);
			fila_iesimaCr(U,L,A,i,n);
		}
		
		printf("Matriz Original: \n\n");
		imp_mat(A,n);
		
		printf("L: \n\n");
		imp_mat(L,n);
		
		printf("U: \n\n");
		imp_mat(U,n);
		
		system("pause");
		
		lib(A,n);
		lib(L,n);
		lib(U,n);
		
	}
}

int jacobi(){
	int dim;
	float norma(float vector1[],float vector2[]);
	float suma_jacobi(float Matriz[], float vector[], int componente);
    int i,j,iteraciones=0;
 
	
    float error,epsilon;
    printf("\n METODO DE JACOBI DE RESOLUCION DE SISTEMAS Ax=b \n");

    printf("Dimension de la matriz A: ");
    scanf("%d",&dim);
    float A[dim][dim],b[dim],x[dim],x_prev[dim],aux[dim];

    printf("\n Elementos de la matriz A: \n");
    for(i=1;i<dim;i++) for(j=1;j<dim;j++){
        printf("A(%d,%d)=",i,j); scanf("%f",&A[i][j]);
    }

    printf("\n Elementos del vector b: \n");
    for(i=1;i<dim;i++){
        printf("b(%d)=",i); scanf("%f",&b[i]);
    }

    printf("\n Error de parada: \n");
    printf("E=",i); scanf("%f",&epsilon);
    error=epsilon+1;
    
    //cominezo algoritmo de Jacobi
    //Error se mide como la norma del vector diferenceia entre la iteracion i e i+1
    printf("\n Valor inicail de la iteracion: \n");
    for(i=1;i<dim;i++){
        printf("x0(%d)=",i); scanf("%f",&x_prev[i]);
    }
    while (error>epsilon){
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++) aux[j]=A[i][j];
            x[i]=(1/A[i][i])*(b[i]-suma_jacobi(aux,x_prev,i));
        }
        error=norma(x,x_prev);

        printf("\n\n Iteracion %d: \n",iteraciones);
        for(i=0;i<dim;i++){
            x_prev[i]=x[i];
            printf("X(%d)=%f \n",i,x[i]);
        }

        iteraciones++;
        if (iteraciones==10) error=epsilon-1;
    }

    printf("Solucion del sistema\n");
    printf("Numero de iteraciones: %d \n", iteraciones);
    for(i=0;i<dim;i++){
        printf("x(%d)=%f\n",i,x[i]);
    }
    return 1;
    system("pause");
}

float norma(float vector1[],float vector2[]){
    float aux=0;
    int i,dim;
    for(i=1;i<dim;i++){
        aux=aux+(vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
    }
    return aux;
}

float suma_jacobi(float Matriz[], float vector[], int componente)
{
    float aux=0;
	int dim;
    int i;
    for(i=0;i<dim;i++){
        if (componente!=i){
            aux=aux+Matriz[i]*vector[i];
        }
    }
    return aux;
}

float ABS(float x){
	if(x < 0)
		return -1*(x);
	else 
		return x; 
}

void sol(float* x,int r){
	int i,n;
	
	printf("%i -> ",r);
	
	for(i = 0; i < n; i++){
		printf("%.6f ",x[i]);
	}
}

void dest(float** A){
	int i,n;
	
	for(i = 0; i < n; i++){
		free(A[i]);
	}
	free(A);
}

void intercambioS(float** A,int r1, int r2){
	float* aux;
	
	aux = A[r1];
	A[r1] = A[r2];
	A[r2] = aux;
	
}

void pedir(float** A){
	int i,j,n;
	
	for(i = 1; i < n; i++){
		printf("Para la ecuacion %i : \n\n",i+1);
		for(j = 1; j < n; j++){
			printf("\tIntroduce a%i : ",j);
			scanf("%f",&A[i][j]);
		}
		system("cls");
	}
	
	printf("Introduccion de terminos independientes\n\n");
	system("pause");
	system("cls");
	
	for(i = 1; i < n;i++){
		printf("Para la ecuacion %i : ",i+1);
		scanf("%f",&A[i][n]);
	}
}

void imprimir(float** A){
	int i,j,n;
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(j == 0){
				printf("%.2fx%i ",A[i][j],j);
			}
			else
				printf("+ %.2fx%i ",A[i][j],j);
		}
		printf("= %.2f",A[i][n]);
		printf("\n");
	}
	printf("\n\n");
	system("pause");
	system("cls");
}

int dom(float** A, int r){
	int i,n;
	int flag = 1;//Suponemos que si es dominante
	
	for(i = 0; i < n; i++){
		if(ABS(A[r][i]) > ABS(A[r][r])){
			flag = 0;
			break;
		}
	}
	return flag;
}

int EDD(float** A){
	int i,j,n;
	int flag = 1;//Suponemos que el SEL es Estrictamente Dominante Diagonalmente
	int band = 1;
	int cont = 0;//Cuenta las veces en las que se intercambian los renglones sin que se arregle lo ED
	int inter = 0;
	
	for(i = 0; i < n && cont < n; i++){
		if(inter){
			band = 1;
			cont = 1;
			for(j = 0; j < n-1 && band; j++){
				intercambio(A,i,j);
				
				if(dom(A,i)){
					i = -1;
					band = 0;
					cont = 0;
				}
				else{
					intercambio(A,i,j);
					cont++;
				}
			}
		}
		else if(i == n-1 && !(dom(A,i))){
			band = 1;
			for(j = 0; j < n-1 && band; j++){
				intercambio(A,i,j);
				
				if(dom(A,i)){
					i = j-1;
					band = 0;
					inter = 1;
				}
				else
					intercambio(A,i,j);
			}
		}
		else{
			if(!(dom(A,i))){
				band = 1;
				for(j = 0; j < n && band; j++){
					intercambio(A,i,j);
					
					if(dom(A,i)){
						i = -1;
						band = 0;
						cont = 0;
					}
					else{
						intercambio(A,i,j);
						cont++;
					}
				}
			}
			else{
				band = 0;
			}
		}
		if(band)
			break;			
	}
	
	if(band || cont >= 3)
		return 0;//No es EDD
	else{
		return 1;
	}
	
}

int Error(float* x,float* ant,int e){
	int i = 0,flag = 1,n;
	float err;
	float v = 1/pow(10,e);
	
	for(i = 0; i < n; i++){
		err = ABS((x[i] - ant[i])/x[i]);
		if(err > v){
			flag = 0; //No cumple la tolerancia aun
			break;
		}
	}
	
	return flag;
}

float xi(float** A,int i,float* x){
	int j,n;
	float sum=0,v;
	float r;
	
	for(j = 0; j < n; j++){
		v = x[j];
		if(j != i)
			sum += A[i][j]*v;
	}
	
	r = (A[i][n]-sum)/A[i][i];
	
	return r;
	
}

float** igualarS(float** A){
	float** B;
	int i,j,n;

	B = (float**)malloc(n*sizeof(float*));

	for(i = 0; i < n; i++){
		B[i] = (float*)malloc(n*sizeof(float));
	}
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			B[i][j] = A[i][j];
		}
	}
	
	return B;
	
}

void sumaS(float **A,float esc,int f, int d){
	float *C;
	int i,n;
	
	C = (float*)malloc(n*sizeof(float));
	
	for(i = 0; i < n; i++){
		C[i] = A[f][i];
		C[i] *= esc;
	}
	for(i = 0; i < n; i++){
		A[d][i] += C[i];
	}
	
	
	free(C);
}

int equalS(float* A, float *B){
	int flag = 1,i,n;
	
	for(i = 0; i < n; i++){
		if(A[i] != B[i]){
			flag = 0;
			break;
		}
	}
	
	return flag;
}

int ceroS(float *A){
	int cont = 0,i,n;
	
	for(i = 0; i < n; i++){
		if(A[i] == 0)
			cont++;
	}
	
	if(cont == n)
		return 1; //Todos los elementos de ese vector son 0
	else
		return 0; //No todos los elementos o ninguno de ese vector son 0
}

float detS(float** A){
	float** B;
	float d = 1,alfa;
	int i,j,n;
	
	B = igualarS(A);
	
	for(j = 0; j < n-1; j++){
		if(ceroS(B[j]) || equalS(B[j],B[j+1])){
			d = 0;
			break;
		}
		else if(B[j][j] == 0){
			d *= -1;
			intercambio(B,j,j+1);
			j--;
		}
		else{
			for(i = j; i < n; i++){
				if(i != j){
					alfa = (-1)*(B[i][j])/B[j][j];
					sumaS(B,alfa,j,i);
				}
			}
			d *= B[j][j];
		}
	}
	
	d*= B[n-1][n-1];
	
	dest(B);
	
	return d;
}

int seidel(){
	float** A;
	float* x_i;//Vector de variables que cambia
	float* x;//Vector de variables finales
	int i,n;
	int e;
	float* ant;
	int k = 2;
	
	printf("Introduce el numero de ecuaciones: "); scanf("%i",&n);
	
	A = (float**)malloc(n*sizeof(float*));
	
	for(i = 0; i < n; i++){
		A[i] = (float*)malloc((n+1)*sizeof(float));
	}
	
	x_i = (float*)malloc(n*sizeof(float));
	
	ant = (float*)malloc(n*sizeof(float));
	
	x = (float*)malloc(n*sizeof(float));

	for(i = 0; i < n; i++){
		x_i[i] = 0;
		x[i] = 0;
		ant[i] = 0;
	}
	
	pedir(A);
	imprimir(A);
	
	
	if(detS(A) == 0){
		printf("El sistema de ecuaciones no tiene solucion unica.\n\n");
		system("pause");
		system("cls");
	}
	else if(!(EDD(A))){
		printf("La convergencia no se garantiza por no ser un sistema EDD\n\n");
		system("pause");
		system("cls");
		
		printf("Introduzca los decimales del error (1^-x) : "); scanf("%i",&e);
		
		printf("x1  x2  x3\n\n");
		
		for(i = 0; i < n; i++){
			x[i] = xi(A,i,x_i);
			x_i[i] = x[i];
		}
		sol(x,1);
		printf("\n");
		do{
			for(i = 0; i < n; i++){
				ant[i] = x[i];
				x[i] = xi(A,i,x_i);
				x_i[i] = x[i];
			}
			sol(x,k);
			k++;
			printf("\n");
		}while(!Error(x,ant,e));
		
		printf("\n\n");
		for(i = 0; i < n; i++){
			printf("x%i = %.6f\n",i,x[i]);
		}
		
		printf("\n\n");
		system("pause");
		system("cls");
	}
	else{
		printf("Introduzca los decimales del error (1^-x) : "); scanf("%i",&e);
		
		printf("x1  x2  x3\n\n");
		
		for(i = 0; i < n; i++){
			x[i] = xi(A,i,x_i);
			x_i[i] = x[i];
		}
		sol(x,1);
		printf("\n");
		do{
			for(i = 0; i < n; i++){
				ant[i] = x[i];
				x[i] = xi(A,i,x_i);
				x_i[i] = x[i];
			}
			sol(x,k);
			k++;
			printf("\n");
		}while(!Error(x,ant,e));
		
		printf("\n\n");
		for(i = 0; i < n; i++){
			printf("x%i = %.6f\n",i,x[i]);
		}
		
		printf("\n\n");
		system("pause");
		system("cls");
		
	}
	
	dest(A);
	free(x);
	free(x_i);
	free(ant);
	
	return 0;
}

typedef float** MD;
typedef float* pentero;

int dimensionMatriz();
MD generaMatriz(int dim);
void liberarMatriz(MD p, int dim);
float cofactor(MD matriz, int orden, int fila, int columna);

MD generaMatriz(int dim){
    MD aux;
    int i;
    if((aux = (MD)malloc(dim*sizeof(pentero))) == NULL){
        return aux;
    }
    
    for( i=0; i<dim; i++){
        if((aux[i] = (pentero)malloc(dim*sizeof(int))) == NULL){
            return aux;
        }
    }
    return aux;
}

void liberarMatriz(MD p, int dim){
	int i;
    for(i=0; i<dim; i++){
        free(p[i]);
    }
    free(p);
}

pentero generaVector(int dim){
    pentero V;
    V = (pentero)malloc(dim*sizeof(float));
    return V;
}

void liberaVector(pentero V){
    free(V);
}

MD inversa(MD mat, int ren){
    int dim, i, j, z, h, k;
    dim = ren;
    float matriz[dim][dim*2], dif, mult;/**Se declara un arreglo bidimension estatico con dimensiones de la matriz
                                        original**/
    for(i=0; i<dim; i++){
        for(j=0; j<dim*2; j++){
            matriz[i][j]=0;/*se inicializa la matriz en "0"*/
        }
    }
    j=dim;
    /*Se transforma la matriz en la indentidad*/
    for(i=0; i<dim; i++){
        matriz[i][j]=1;/*1´s en la diagonal principal*/
        j++;
    }
    /*Se copia la matriz original en la matriz estatica*/
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            matriz[i][j]=mat[i][j];
        }
    }
    /*Gaus-Jordan*/
    printf("\n");
    for(z=0; z<dim; z++){
        dif = matriz[z][z];
        for(h=0; h<dim*2; h++){
            matriz[z][h]=matriz[z][h]/dif;/*Se divide cada fila entre pivote correspondiente */
        }
        for(i=0; i<dim; i++){
            if (i!=z){
                mult= -matriz[i][z];
                for (j=0; j<dim*2; j++){
                    matriz[i][j]=matriz[z][j]*mult + matriz[i][j];/*"matriz" se transforma en la indentidad*/
                }
            }
        }
    }
    /*Se copia la "matriz" en "mat" ya que esta se ransforma en la inversa de "mat"*/
    for(i=0; i<dim; i++){
        k=ren;
        for (j= 0; j<dim; j++){
            mat[i][j]= matriz[i][k];
            k++;
        }
    }
    return mat;
}

/**Determinante de la matriz*/
float determinante(MD matriz, int orden){
    float det;
    int i;
    if(orden == 1){
        det = matriz[0][0];
    }else{
        for(i=0; i<orden; i++){
            det = det + matriz[0][i]*cofactor(matriz, orden, 0, i);/*Se llama la funcion "cofactor" mientras
                                                                       la dimension de "matriz" sea distinto de 1**/
        }
    }
    return det;/*Se retorna el determinante cuando la dimension de "matriz" sea 1*/
}

/*****Cofactor***/
float cofactor(MD matriz, int orden, int fila, int columna){
    MD subB;
    subB = generaMatriz(orden);
    int n = orden-1;
    int x = 0;
    int y = 0;
    int i,j;
    /*recorrido interno*/
    for(i=0; i<orden; i++){
        for(j=0; j<orden; j++){
            if(i!=fila && j!=columna){
                subB[x][y]=matriz[i][j];
                y++;
                if(y >= n){
                    x = x + 1;
                    y = 0;
                }
            }
        }
    }
    return pow(-1.0, fila+columna)*determinante(subB, n);/*Antes de retornar el valor, éste llama a la funcion "determinante"
                                                         hasta que "matriz" sea de dimension 1**/
}

void leerMatriz(MD p, int dim){
    int i,j;
	for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf("\t p[%d][%d]: ",i+1,j+1);
            scanf("%f",&p[i][j]);
        }
    }
}

void mostrarMatriz(MD p, int dim){
    int i, j;
	for(i=0; i<dim; i++){
        printf("\n");
        for(j=0; j<dim; j++){
            printf("\t %f ",p[i][j]);
        }
    }
}

void ingresaVector(pentero V, int dim){
    int i;
	for(i=0; i<dim; i++){
        printf("\n [%i] : ",i+1);
        scanf("%f",&V[i]);
    }
}

void imprimeVector(pentero V, int dim){
    int i;
	for(i=0; i<dim; i++){
        printf("\n %f ",V[i]);
    }
}

void imprimeInv(MD inv, int dim){
    int i,j;
	printf("\n");
    for( i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf(" %f ",inv[i][j]);
        }
        printf("\n\n");
    }
}

void producto(MD mat, pentero vectCons,int dim){
    int i,j,k;
	float sol[dim];
    for(j=0; j<dim; j++){
        float acumula=0.0;
        for(k=0; k<dim; k++){
            acumula+= mat[j][k]*vectCons[k];
        }
        sol[j] = acumula;
    }
    for(i=0; i<dim; i++){
        vectCons[i] = sol[i];
    }
}

int Inver(){
    MD A,B;
    int i,j;
    pentero vectCons;
    int dim=0;
    printf("\n Ingresa los datos de la matriz A ");
    printf("\n Ingresa el rango de la matriz: ");
    scanf("%d",&dim);
    A = generaMatriz(dim);
    B = generaMatriz(dim);
    printf("\n Ingresa los elementos de la matriz \n");
    leerMatriz(A, dim);
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            B[i][j] = A[i][j];
        }
    }
    printf("\n Matriz A \n");
    mostrarMatriz(A, dim);
    int det = determinante(A, dim);
    if(det == 0){
        printf("\n El sistema no tiene solucion. ");
        return 13;
    }
    printf("\n Ingresa el vector de terminos independientes \n");
    vectCons = generaVector(dim);
    ingresaVector(vectCons, dim);
    printf("\n\n Matriz inversa de A");
    B = inversa(B, dim);
    imprimeInv(B, dim);
    printf("\n Vector de terminos independientes ");
    imprimeVector(vectCons, dim);
    printf("\n\n Vector solucion");
    producto(B, vectCons, dim);
    imprimeVector(vectCons, dim);
    printf("\n\n Comprobacion: A*Vector solucion = Vector de terminos independientes ");
    producto(A, vectCons, dim);
    imprimeVector(vectCons, dim);
    liberarMatriz(A, dim);
    liberaVector(vectCons);
    liberarMatriz(B, dim);
}

void imp_zero(int** A, int n){
	int i,j;
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("[%i] ",A[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}

void col_inter(float** B, int c1, int c2, int n){
	float aux = 0;
	int i,j;
	
	for(i = 0; i < n; i++){
		aux = B[i][c1];
		B[i][c1] = B[i][c2];
		B[i][c2] = aux;
	}
	
}

void col_pivote(float** A, int r, int c, int n){
	float p = A[r][c];
	int i;
	
	for(i = 0; i < n; i++){
		if(i != r)
			A[i][c] /= p;
	}
	
}

void ren_pivote(float** A, int r, int c, int n){
	float p = A[r][c];
	int i;

	for(i = 0; i < n; i++){
		if(i != c)
			A[r][i] /= (-1*p);
	}
	
}

void vecinos(float** A, int r, int c, int n){
	int i,j;
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i != r && j != c){
				A[i][j] += A[r][j]*A[i][c];
			}
		}
	}
}

void intercambioIn(float** A, int r, int c, int n){
	ren_pivote(A,r,c,n);
	vecinos(A,r,c,n);
	col_pivote(A,r,c,n);
	A[r][c] = 1/A[r][c];
}

float detI(float** A,int n){//Funcion que retorna el determinante de una matriz
	float **B;
	float d = 1;
	float alfa;
	int j,i;
	
	B = igualar(A,n);
	
	for(j = 0; j < n-1; j++){
		if(cero(B[j],n) || equal(B[j],B[j+1],n)){
			d = 0;
			break;
		}
		else if(B[j][j] == 0){
			d *= -1;
			intercambioIn(B,j,j+1,n);
			j--;
		}
		else{
			for(i = j; i < n; i++){
				if(i != j){
					alfa = (-1)*(B[i][j])/B[j][j];
					suma(B,alfa,j,i,n);
				}
			}
			d *= B[j][j];
		}
	}
	d*= B[n-1][n-1];
	lib(B,n);
	
	return d;
}

void uno(int** zero, int r, int c, int n){
	int i,j;
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == r || j == c)
				zero[i][j] = 1;
		}
	}
}

int mayor(float* A, int* zero, int n){
	int i;
	int c = -1;
	float max = 0;
	
	for(i = 0; i < n; i++){
		if(!zero[i]){
			if(i == n-1)
				c = i;
			else{
				if(ABS(A[i]) > ABS(max)){
					max = A[i];
					c = i;
				}
			}
		}
	}
	
	return c;
}

void imp_vect(float* a, int n){
	int i;
	
	for(i = 0; i < n; i++){
		if(i == 0){
			if(a[0] != 0)
				printf("%.3f x%i ",a[i],i);
			else
				printf("\t");
		}
		else{
			if(a[i] < 0)
				printf("- %.3f ",ABS(a[i]));
			else if(a[i] > 0)
				printf("+ %.3f ",a[i]);
			else
				printf("\t");
				
			if(a[i] != 0)
				printf("x%i ",i);
		}
	}
}

void imp_SEL(float** A, float* b, int n){
	int i,j;
	
	for(i = 0; i < n; i++){
		imp_vect(A[i],n);
		printf("= %.3f\n",b[i]);
	}
	printf("\n\n");
}

float* sol_sistema(float** A, float* b, int n){
	int i,j;
	float* x;
	float sum = 0;
	
	do{
		x = (float*)calloc(n,sizeof(float));
	}while(x == NULL);
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			sum += A[i][j]*b[j];
		}
		x[i] = sum;
		sum = 0;
	}
	
	return x;
}

int Int(){
	float** A;
	float** B;
	float* sol;
	int** zero;
	float* b;
	int* x;
	int* y;
	int n;
	int i,j;
	int pivote;
	int aux;
	
	printf("Introduce el numero de incognitas del sistema : "); scanf("%i",&n);
	system("cls");
	
	printf("Introduce los valores respectivos de la matriz de coeficientes: \n\n");
	
	A = creaMatriz(n);
	ingresaMatriz(A,n);
	
	if(!detI(A,n)){
		printf("El determinante del sistema es cero, por lo que el sistema no es compatible determinado.\n\n");
		system("pause");
		system("cls");
	}
	else{
		B = igualar(A,n);
		
		do{
			b = (float*)calloc(n,sizeof(float));
		}while(b == NULL);
		
		system("cls");
		printf("Ahora introduce los valores del vector independiente.\n\n");
		for(i = 0; i < n; i++){
			printf("b%i = ",i+1);
			scanf("%f",&b[i]);
		}
		
		printf("\n\n");
		system("pause");
		system("cls");
		
		printf("La matriz original es: \n\n");
		
		imp_mat(A,n);
		
		do{
			zero = (int**)malloc(n*sizeof(int*));
		}while(zero == NULL);
		
		for(i = 0; i < n; i++){
			do{
				zero[i] = (int*)calloc(n,sizeof(int));
			}while(zero[i] == NULL);
		}
		do{
			x = (int*)calloc(n,sizeof(int));
		}while(x == NULL);
		
		for(i = 1; i < n; i++){
			x[i] = i;
		}
		
		do{
			y = (int*)calloc(n,sizeof(int));
		}while(x == NULL);
		
		for(i = 1; i < n; i++){
			y[i] = i;
		}
		printf("Aplicando el metodo del intercambio: \n\n");
		
		for(i = 0; i < n; i++){
			printf("Iteracion %i\n\n",i+1);
			pivote = mayor(A[i],zero[i],n);
			x[i] = pivote;
			y[pivote] = i;
			intercambioIn(A,i,pivote,n);
			uno(zero,i,pivote,n);
			imp_mat(A,n);
		}
		
		system("pause");
		system("cls");
		
		printf("Ordenando la matriz inversa de coeficientes con el orden normal:\n\n");
		
		for(i = 0; i < n; i++){
			for(j = 0; j < n-1; j++){
				if(x[j] > x[j+1]){
					aux = x[j];
					x[j] = x[j+1];
					x[j+1] = aux;
					intercambioIn(A,j,j+1,n);
				}
				if(y[j] > y[j+1]){
					aux = y[j];
					y[j] = y[j+1];
					y[j+1] = aux;
					col_inter(A,j,j+1,n);
				}
			}
		}
		
		imp_mat(A,n);
		
		system("pause");
		system("cls");
		
		sol = sol_sistema(A,b,n);
		
		printf("Para el sistema de ecuaciones: \n\n");
		imp_SEL(B,b,n);
		
		printf("La solucion es: \n\n");
		for(i = 0; i < n; i++){
			printf("\tx%i = %.3f\n",i+1,sol[i]);
		}
		
		printf("\n\n");
		system("pause");
		system("cls");
		
		for(i = 0; i < n; i++){
			free(zero[i]);
		}
		free(zero);
		free(x);
		free(y);
		lib(B,n);
	}
	
	lib(A,n);
	return 0;
	
	
}

int bis(void){
	
	
    int iter,xmax,res;
   float x,tol,fxa,fxb,test;
   float xa,xb,xn,xnold,ea,xa1,xb1;
	printf("Para la ecuacion: Sen(x)Cos(x)");
	do{
	   printf("\nIntroduce el valor \''a\'' de tu intervalo: ");
   scanf("%f",&xa);

      printf("\nIntroduce el valor \''b\'' de tu intervalo: ");
   scanf("%f",&xb);

    printf("\nIngrese el numero maximo de interaciones:  ");
   scanf("%d",&xmax);

    printf("Ingrese la tolerancia: ");
   scanf("%f",&tol);

   iter=0;
    xn=0;
		system("cls");
		res=RAIZ(xa)*RAIZ(xb);
		printf("%i",res);
		//printf("\n\n\nEl metodo no converge a una raiz real.\n");
        //printf("introduzca otro intervalo");
		}while(res>0);
		
 	 if(RAIZ(xa)*RAIZ(xb)<0){

    printf("\nX          [a, b]              Xn         Error \n");

    do{

        iter++;


        fxa= RAIZ(xa);

        fxb=RAIZ(xb);

        xb1=xb;

        xa1=xa;

        xnold=xn;

        xn=(xa+xb)/2;

        ea=fabs((xn-xnold)/xn);

        test=RAIZ(xn)*RAIZ(xa);

        if(test<0)
            xb=xn;

        else if(test >0)
            xa=xn;

        else
            ea=0;

        printf("\n%d  [%lf, %lf]    %lf    %lf ",iter,xa1,xb1,xn,ea);


    } while(ea>tol && iter<xmax );

}
    printf("\n\n\tLa raiz aprox de su funcion es:  %lf \n\n\n", xn);


   return 0;
}

int Fp(void){
    int iter,xmax,res;
   float x,tol,fxa,fxb,test;
   float xa,xb,xn,xnold,ea,xa1,xb1;
   printf("Para la ecuacion:    (e^x)Cos(x)) ");
do{
   printf("\nIntroduce el valor \''a\'' de tu intervalo: ");
   scanf("%f",&xa);

      printf("\nIntroduce el valor \''b\'' de tu intervalo: ");
   scanf("%f",&xb);

    printf("\nIngrese el numero maximo de interaciones:  ");
   scanf("%d",&xmax);

    printf("Ingrese la tolerancia: ");
   scanf("%f",&tol);
	 iter=0;
    xn=0;
	
		system("cls");
		res=RAIZ1(xa)*RAIZ1(xb);
      //  printf("\n\n\nEl metodo no converge a una raiz real.\n");
       // printf("escoja otro intervalo otro intervalo");
  
    }while(res>0);

  if(RAIZ1(xa)*RAIZ1(xb)<0){

    printf("\nX          [a, b]              Xn         Error \n");

    do{

        iter++;


        fxa= RAIZ1(xa);

        fxb=RAIZ1(xb);

        xb1=xb;

        xa1=xa;

        xnold=xn;

        xn=(xa*(RAIZ1(xb))-xb*(RAIZ1(xa)))/(RAIZ1(xb)-RAIZ1(xa));

        ea=fabs((xn-xnold)/xn);

        test=RAIZ1(xn)*RAIZ1(xa);

        if(test<0)
            xb=xn;

        else if(test >0)
            xa=xn;

        else
            ea=0;

        printf("\n%d  [%lf, %lf]    %lf    %lf ",iter,xa1,xb1,xn,ea);


    } while(ea>tol && iter<xmax );

}

    printf("\n\n\tLa raiz aprox de su funcion es:  %lf \n\n\n",xn);


   return 0;
}


typedef float** MATRIZ;
typedef float* vectorConst;

MATRIZ generaMatrizR(int dim);
void liberaMatrizR(MATRIZ A, int dim);
vectorConst generaVectorR(int dim);
void liberaVectorR(vectorConst V);
void ingresaMatrizR(MATRIZ A, int dim);
void ingresaVectorR(vectorConst V, int dim);
void imprimeMatrizR(MATRIZ A, int dim);
void imprimeVectorR(vectorConst V, int dim);
void transformaSistemaR(MATRIZ A, MATRIZ B, vectorConst V, vectorConst V2, int dim);
void obtieneRaizR(MATRIZ B, vectorConst V2, vectorConst vect, int dim);


MATRIZ generaMatrizR(int dim){
	 int i;
    MATRIZ aux;
    if((aux = (MATRIZ)malloc(dim*sizeof(pentero))) == NULL){
        return 0;
    }  
   
	for(i=0; i<dim; i++){
        if((aux[i] = (pentero)malloc(dim*sizeof(float))) == NULL){
            return 0;
        }
    }
    return aux;
}

void liberaMatrizR(MATRIZ A, int dim){
    int i;
	for(i=0; i<dim; i++){
        free(A[i]);
    }
    free(A);
}
vectorConst generaVectorR(int dim){
    vectorConst V;
    V = (vectorConst)malloc(dim*sizeof(float));
    return V;
}
void liberaVectorR(vectorConst V){
    free(V);
}
void ingresaMatrizR(MATRIZ A, int dim){
	int i,j;
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf("\n [%d][%d] : ",i+1,j+1);
            scanf("%f",&A[i][j]);
        }
    }
}
void ingresaVectorR(vectorConst V, int dim){
    int i;
	for(i=0; i<dim; i++){
        printf("\n [%i] : ",i+1);
        scanf("%f",&V[i]);
    }
}
void imprimeMatrizR(MATRIZ A, int dim){
    int i,j;
	printf("\n");
    for(i=0; i<dim; i++){
        for(j=0; j<dim ; j++){
            printf(" %.4f ",A[i][j]);
        }
        printf("\n");
    }
}
void imprimeVectorR(vectorConst V, int dim){
    int i;
	for(i=0; i<dim; i++){
        printf("\n %.4f ",V[i]);
    }
}
void transformaSistemaR(MATRIZ A, MATRIZ B, vectorConst V, vectorConst V2, int dim){
    int i,j;
	for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            B[i][j] = (-1)*A[i][j]/A[i][i];
            if(j == dim-1){
                V2[i] = (-1)*V[i]/A[i][i];
            }
        }
    }
}
float productR(MATRIZ B, vectorConst V2, vectorConst vect, int dim, int i);

void obtieneRaizR(MATRIZ B, vectorConst V2, vectorConst vect, int dim){
    float raiz;
    float arreglo[dim];
    int itera=0;
    float error = 0.00001;
    int i;
	for(i=0; i<dim; i++){
        vect[i] = 0;
        arreglo[i] = 0;
    }
    int indi = 0;
    int and;
    printf("\n| Iteracion|");
    printf(" ");
    for(i=0; i<dim; i++){
        printf("   X%d    |",i+1);
    }
    printf("      |");
    for(i=0; i<dim; i++){
        printf("     R%d    |",i+1);
    }
    printf("       |  MAXIMO |");
    do{
        printf("\n\n|    %d     | ",indi);
        indi+=1;
        for(i=0; i<dim; i++){
            printf(" %.5f |",vect[i]);
        }
        printf("      |");
        float prod=0;
        for(i=0; i<dim; i++){
            arreglo[i] = productR(B, V2, vect, dim, i);
        }
        for(i=0; i<dim; i++){
            printf("  %.5f  |",arreglo[i]);
        }
        float maxi = 0.0;
        for(i=0; i<dim; i++){
            if(fabs(arreglo[i]) > maxi){
                maxi = fabs(arreglo[i]);
            }
        }
        printf("       | %.5f |",maxi);
        for(i=0; i<dim; i++){
            if(fabs(arreglo[i]) == maxi){
                vect[i] = vect[i]+arreglo[i];
            }else{
                vect[i] = vect[i];
            }
        }
        and=0;
        int i;
        for(i=0; i<dim; i++){
            if(fabs(arreglo[i]) < error){
                and+=1;
            }
        }
    }while(and != dim);
    printf("\n Se detiene el metodo...");
}
float productR(MATRIZ B, vectorConst V2, vectorConst vect, int dim, int i){
	int j;
    float prod=0.0;
    for(j=0; j<dim; j++){
        prod += B[i][j]*vect[j];
    }
    prod = prod-V2[i];
    return prod;
}

int rel(){
    int dim;
    float raiz;
    MATRIZ A, B;
    vectorConst V, V2, vect;
    printf("\n Ingresa el rango de la matriz: ");
    scanf("%d",&dim);
    A = generaMatrizR(dim);
    B = generaMatrizR(dim);
    V = generaVectorR(dim);
    V2 = generaVectorR(dim);
    vect = generaVectorR(dim);
    printf("\n Ingresa los datos de la matriz ");
    ingresaMatrizR(A, dim);
    printf("\n Ingresa los datos del vector: ");
    ingresaVectorR(V, dim);
    printf("\n Matriz ingresada ");
    imprimeMatrizR(A, dim);
    printf("\n Vector ingresado ");
    imprimeVectorR(V, dim);
    transformaSistemaR(A, B, V, V2, dim);
    printf("\n Matriz transformada. ");
    imprimeMatrizR(B, dim);
    printf("\n Vector transformado. ");
    imprimeVectorR(V2, dim);
    obtieneRaizR(B, V2, vect, dim);


   //Se libera memoria de las matrices y vectores.
    liberaMatrizR(A, dim);
    liberaMatrizR(B, dim);
    liberaVectorR(V);
    liberaVectorR(V2);
    liberaVectorR(vect);
}

int sec(void){
    int iter,xmax,res;
   float x,tol,fxa,fxb,test;
   float xa,xb,xn,xnold,ea;
   printf("Para la ecuacion: (x-5)(x+3)Sen^3 (x)");

	do{
   printf("\nIntroduce el valor inicial: ");
   scanf("%f",&xa);

      printf("\nIntroduce el valor final: ");
   scanf("%f",&xb);


    printf("Ingrese la tolerancia: ");
   scanf("%f",&tol);

   iter=0;
    xn=0;
	system("cls");
    res=RAIZ3(xa)*RAIZ3(xb);
        //printf("\n\n\nEl metodo no converge a una raiz real.\n");
       // return 0;
    }while(res>0);

     if(RAIZ3(xa)*RAIZ3(xb)<0){

    printf("\nI    Xn-1        Xn       Xn+1       Error \n",tol);


    do{

        iter++;


        fxa= RAIZ3(xa);

        fxb=RAIZ3(xb);

        xnold=xn;

        xn=xb-(fxb*(xa-xb))/(fxa-fxb);

        ea=fabs((xn-xnold)/xn);

        test=fabs(RAIZ3(xn));




        printf("\n%d  %lf   %lf   %lf    %lf ",iter,xa,xb,xn,ea);

        xa=xb;
        xb=xn;


    } while(test>tol );

}

    else{

        printf("Error, vuelva a intentarlo");
    }

    printf("\n\n\tLa raiz aprox de su funcion es:  %lf \n\n\n", xn);


   return 0;
}

int newt(void){
    int iter,xmax;
   float x,tol,fxa,fxb,test;
   float xa,xn,xnold,ea;
printf("Para la ecuacion:   (x^2)Cos^2 (x) ");

   printf("\nIntroduce el valor inicial: ");
   scanf("%f",&xa);


    printf("Ingrese la tolerancia: ");
   scanf("%f",&tol);

   iter=0;
    xn=0;


    printf("\nI    Xn-1        Xn       Xn+1       Error \n",tol);


    do{

        iter++;


        fxa= RAIZ2(xa);

        fxb=DER(xa);

        xnold=xn;

        xn=xa-(RAIZ2(xa))/(DER(xa));

        ea=fabs((xn-xnold)/xn);

        test=fabs(RAIZ2(xn));




        printf("\n%d  %lf   %lf   %lf    %lf ",iter,fxa,fxb,xn,ea);

        xa=xn;


    } while(test>tol );



    printf("\n\n\tLa raiz aprox de su funcion es:  %lf \n\n\n", xn);


   return 0;
}

int dimensionMatrizJ();
MD generaMatrizJ(int dim);
MD generaMatrizDdimJ(int fila, int col);
void liberarMatrizJ(MD p, int dim);
void liberarMatrizDdimJ(MD p, int fila, int col);
float cofactorJ(MD matriz, int orden, int fila, int columna);
/*********** Se genera memoria para una matriz cuadrada **************/
MD generaMatrizJ(int dim){
	int i;
    MD aux;
    if((aux = (MD)malloc(dim*sizeof(pentero))) == NULL){
        return 0;
    }
    for(i=0; i<dim; i++){
        if((aux[i] = (pentero)malloc(dim*sizeof(int))) == NULL){
            return 0;
        }
    }
    return aux;
}
/******************Se genera memoria del monticulo para una matriz no cuadrada *********************/
MD generaMatrizDdimJ(int fila, int col){
	int i;
    MD aux;
    if((aux = (MD)malloc(fila*sizeof(pentero))) == NULL){
        return 0;
    }
    for(i=0; i<fila; i++){
        if((aux[i] = (pentero)malloc(col*sizeof(int))) == NULL){
            return 0;
        }
    }
    return aux;
}
/**************** Se libera memoria ******************/
void liberarMatrizJ(MD p, int dim){
    int i;
	for(i=0; i<dim; i++){
        free(p[i]);
    }
    free(p);
}

void liberarMatrizDdimJ(MD p, int fila, int col){
    int i;
	for( i=0; i<fila; i++){
        free(p[i]);
    }
    free(p);
}
/****************  Se genera un arreglo unidimensional *******************/
pentero generaVectorJ(int dim){
    pentero V;
    V = (pentero)malloc(dim*sizeof(float));
    return V;
}
/*************** Liberacion de memoria del arreglo *******************/
void liberaVectorJ(pentero V){
    free(V);
}
/****************** Se calcula la inversa de una matriz ************************/
MD inversaJ(MD mat, int ren){
    int dim, i, j, z, h, k;
    dim = ren;
    float matriz[dim][dim*2], dif, mult;//Se declara un arreglo bidimension estatico con dimensiones de la matriz
                                        //original**/
    for(i=0; i<dim; i++){
        for(j=0; j<dim*2; j++){
            matriz[i][j]=0;//se inicializa la matriz en "0"**/
        }
    }
    j=dim;
    //Se transforma la matriz en la indentidad**/
    for(i=0; i<dim; i++){
        matriz[i][j]=1;//1´s en la diagonal principal**/
        j++;
    }
    //Se copia la matriz original en la matriz estatica**/
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            matriz[i][j]=mat[i][j];
        }
    }
    //Gaus-Jordan
    printf("\n");
    for(z=0; z<dim; z++){
        dif = matriz[z][z];
        for(h=0; h<dim*2; h++){
            matriz[z][h]=matriz[z][h]/dif;//Se divide cada fila entre pivote correspondiente **/
        }
        for(i=0; i<dim; i++){
            if (i!=z){
                mult= -matriz[i][z];
                for (j=0; j<dim*2; j++){
                    matriz[i][j]=matriz[z][j]*mult + matriz[i][j];//"matriz" se transforma en la indentidad**/
                }
            }
        }
    }
    //Se copia la "matriz" en "mat" ya que esta se ransforma en la inversa de "mat"**/
    for(i=0; i<dim; i++){
        k=ren;
        for (j= 0; j<dim; j++){
            mat[i][j]= matriz[i][k];
            k++;
        }
    }
    return mat;
}
/***Determinante de la matriz**/
float determinanteJ(MD matriz, int orden){
	int i;
    float det;
    if(orden == 1){
        det = matriz[0][0];
    }else{
        for(i=0; i<orden; i++){
            det = det + matriz[0][i]*cofactor(matriz, orden, 0, i);//Se llama la funcion "cofactor" mientras
                                                                   //la dimension de "matriz" sea distinto de 1**/
        }
    }
    return det;//Se retorna el determinante cuando la dimension de "matriz" sea 1**/
}
/************ Cofactor **********/
float cofactorJ(MD matriz, int orden, int fila, int columna){
    MD subB;
    subB = generaMatriz(orden);
    int n = orden-1;
    int x = 0;
    int y = 0;
    int i,j;
    //recorrido interno
    for(i=0; i<orden; i++){
        for(j=0; j<orden; j++){
            if(i!=fila && j!=columna){
                subB[x][y]=matriz[i][j];
                y++;
                if(y >= n){
                    x = x + 1;
                    y = 0;
                }
            }
        }
    }
    return pow(-1.0, fila+columna)*determinante(subB, n);//Antes de retornar el valor, éste llama a la funcion "determinante"
                                                         //hasta que "matriz" sea de dimension 1**/
}
void leerMatrizJ(MD p, int dim){
    int i,j;
	for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf("\t p[%d][%d]: ",i+1,j+1);
            scanf("%f",&p[i][j]);
        }
        printf("\n");
    }
}
void mostrarMatrizJ(MD p, int dim){
    int i,j;
	for(i=0; i<dim; i++){
        printf("\n");
        for(j=0; j<dim; j++){
            printf("\t %f ",p[i][j]);
        }
    }
}
void ingresaVectorJ(pentero V, int dim){
    int i; 
	for(i=0; i<dim; i++){
        printf("\n [%i] : ",i+1);
        scanf("%f",&V[i]);
    }
}
void imprimeVectorJ(pentero V, int dim){
    int i;
	for( i=0; i<dim; i++){
        printf("\n %f ",V[i]);
    }
}
void imprimeInvJ(MD inv, int dim){
    int i,j;
	printf("\n");
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf(" %f ",inv[i][j]);
        }
        printf("\n\n");
    }
}
/************* Se imprime la matriz no cuadrada *************/
void mostrarMatrizDdimJ(MD mat, int fila, int col){
    int i,j;
	for(i=0; i<fila; i++){
        printf("\n");
        for(j=0; j<col; j++){
            printf(" %f ",mat[i][j]);
        }
    }
}
/************** Multiplicacion por inversas ******************/
void multiInverJ(MD mat1, MD mat2, int ren2, int col2, int ren1, int col1){
    int fila=ren1, colu=col2;
    int i,j,k;
	float aux[fila][colu]; // Se declara un arreglo bidimensional en el cual se almacenara el resultado final ******/
    float acumula;
    for(i=0; i<ren1; i++){
        for( j=0; j<col2; j++){
            acumula=0;
            for(k=0; k<col1; k++){
                acumula += mat2[i][k]*mat1[k][j];
            }
            aux[i][j] = acumula;
        }
    }
    for( i=0; i<fila; i++){
        for( j=0; j<colu; j++){
            mat1[i][j] = aux[i][j];
        }
    }
}
/************* Producto de dos matrices no cuadradas ******************************/
void multiJ(MD matCopy, MD matMult, int ren1, int col1, int ren2, int col2){
    int fila=ren1, colu=col2;
    float aux[fila][colu]; // Se declara un arreglo bidimensional en el cual se almacenara el resultado final ******/
    float acumula;
    int i,j,k;
    // Se realiza el producto de matrices no cuadradas, solo con (columnas de A) = (filas de B) ********/
    for(i=0; i<ren1; i++){
        for(j=0; j<col2; j++){
            acumula=0;
            for(k=0; k<col1; k++){
                acumula += matCopy[i][k]*matMult[k][j];
            }
            aux[i][j] = acumula;
        }
    }
    // Se almacena el producto en mat1 ****************/
//    printf("\n Se imprime desde multiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii \n");
    for(i=0; i<fila; i++){
        for(j=0; j<colu; j++){
            matCopy[i][j] = aux[i][j];
//            printf(" %f ",matCopy[i][j]);
        }
//        printf("\n");
    }
}
/************* Producto entre una matriz y un vector unidimensional *****************************/
void productoMatVecJ(MD mat, pentero vector, int fila, int col){
    float aux[fila];
    float acum=0;
    int i,j;
    for(i=0; i<fila; i++){
        for(j=0; j<col; j++){
            acum += mat[i][j]*vector[j];
        }
        aux[i] = acum;
        acum = 0;
    }
//    printf("\n Se imprime desde vectorrrrrrrrrrrrrrrrrrrrrrrrr \n");
    for(i=0; i<fila; i++){
        vector[i] = aux[i];
//        printf("\n %f ",vector[i]);
    }
}
/**************** Copia matriz ********************/
void copiaMatrizJ(MD mat1, MD mat2, int fila, int col){
    int i,j;
	for(i=0; i<fila; i++){
        for(j=0; j<col; j++){
            mat1[i][j] = mat2[i][j];
        }
    }
}
/**************** Suma de matrices no cuadradas ******************/
void sumaProductJ(MD matFin, MD matCopy, MD matMult,int ren1, int col1, int ren2, int col2, int ren3, int col3){
    MD matCopy1;
    int i,j;
    matCopy1 = generaMatrizDdimJ(ren2, col2);
    copiaMatrizJ(matCopy1, matCopy, ren2, col2);
    multiJ(matCopy, matMult, ren2, col2, ren3, col3);
    for(i=0; i<ren1; i++){
        for(j=0; j<col1; j++){
            matFin[i][j] -= matCopy[i][j];
        }
    }
    copiaMatrizJ(matCopy, matCopy1, ren2, col2);
    liberarMatrizDdimJ(matCopy1, ren2, col2);
}
/******************* Se realiza la siguiente operacio: vector 2 = -(A21 * vector 1)+vector 2 ***************/
void sumaProducVectorJ(pentero vectorFin, MD copy, pentero vectorPro, int ren, int col){
    pentero vectorProCopy; //Se crea un arreglo unidimensional donde se almacenara el arreglo "copy" para que éste
       int i;                     //no pierda sus valores.
    vectorProCopy = generaVector(col);
    for(i=0; i<col; i++){
        vectorProCopy[i] = vectorPro[i];
//        printf("\n %f ",vectorProCopy[i]);
    }
//    printf("\n \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----------\n");
    productoMatVecJ(copy, vectorProCopy, ren, col);
    for(i=0; i<ren; i++){
        vectorFin[i] -= vectorProCopy[i];
//        printf("\n %f ",vectorFin[i]);
    }
    liberaVector(vectorProCopy);//Se libara memoria del arreglo creado en ésta funcion.
}
/*********** Iteracion de Gauss-Jordan particionado *************/
void gausJordanParticionadoJ(MD A,pentero vector,int dim){
    int p=0;
    int i,j;
    int k=0;
    do{
        printf("\n Ingresa el valor de la particion: ");
        scanf("%d",&p);
    }while(p<2 || p>(dim-2));
    //Se genera memoria del monticulo
    MD AuuInv, AddInv, A12, A11, A21, A21copy, A22, A12copy;
    pentero B1, B2, Solucion;
    AuuInv = generaMatriz(p);
    AddInv = generaMatriz(dim-p);
    A11 = generaMatriz(p);
    A12 = generaMatrizDdimJ(p, dim-p);
    B1 = generaVector(p);
    A21 = generaMatrizDdimJ(dim-p, p);
    A21copy = generaMatrizDdimJ(dim-p, p);
    A22 = generaMatrizJ(dim-p);
    B2 = generaVectorJ(dim-p);
    A12copy = generaMatrizDdimJ(p, dim-p);
    Solucion = generaVectorJ(dim);
    //Se almacena en matriz A11 los elementos correspondientes
    for(i=0; i<p; i++){
        for(j=0; j<p; j++){
            A11[i][j] = A[i][j];
            AuuInv[i][j] = A[i][j];
        }
    }
    //Se almacena en matriz A12 los elementos correspondientes
    for( i=0; i<p; i++){
        for(j=0; j<(dim-p); j++){
            A12[i][j] = A[i][j+p];
        }
    }
    //Se almacena en matriz A21 los elementos correspondientes
    for(i=0; i<dim-p; i++){
        k = p;
        for(j=0; j<p; j++){
            A21[i][j] = A[k+i][j];
        }
    }
    //Se almacena en matriz A22 los elementos correspondientes
    for(i=0; i<dim-p; i++){
        for(j=0; j<dim-p; j++){
            A22[i][j] = A[p+i][p+j];
        }
    }
    //Se almacena en vector B1 los elementos correspondientes
    for(i=0; i<p; i++){
        B1[i] = vector[i];
    }
    //Se almacena en vector B2 los elementos correspondientes
    for(i=0; i<dim-p; i++){
        B2[i] = vector[p+i];
    }
    //Se imprime la matriz con la particion correspondiente
    printf("\n ITERACION 0 \n");
    for(i=0; i<p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A11[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A12[i][n]);
                n+=1;
            }else{
                printf("\t  %f  ||",B1[i]);
            }
        }
    }
    printf("\n");
    for(i=0; i<dim-p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A21[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A22[i][n]);
                n+=1;
            }else {
                printf("\t  %f  ||",B2[i]);
            }
        }
    }
    //Termina proceso de almacenamiento
    inversa(AuuInv, p); //Se calcula la inversa de A uno x uno
    multiInverJ(A11, AuuInv, p, p, p, p); //Se realiza la siguiente operacion A11 = A11^1 x A11
    multiInverJ(A12, AuuInv, p, dim-p, p, p); //Se realiza la siguiente operacion: A12 = A11^1 x A12
    copiaMatrizJ(A12copy, A12, p, dim-p);//Se copia A12 * A11^-1 en A12copy
    productoMatVecJ(AuuInv, B1, p, p); //Se realiza la siguiente operacion: B1 =  A11^1 x B1

    //se hace cero A21
    copiaMatrizJ(A21copy, A21, dim-p, p);//Se copia la matriz A21 en A21copy para: A21 =( -1)A21copy * A11 +A21
    sumaProductJ(A21, A21copy, A11, dim-p, p, dim-p, p, p, p);

    //Se hace la operacion A22 = -(A21 * A12) + A22
    sumaProductJ(A22, A21copy, A12, dim-p, dim-p, dim-p, p, p, dim-p);

    //Se realiza la siguiente operacion: B2 = -(A21 * B1) + B2
    sumaProducVectorJ(B2, A21copy, B1, dim-p, p);
    /************************* SIGUIENTE ITERACION **************************************/
    printf("\n ITERACION 1 \n");
    for(i=0; i<p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A11[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A12[i][n]);
                n+=1;
            }else{
                printf("\t  %f  ||",B1[i]);
            }
        }
    }
    printf("\n");
    for(i=0; i<dim-p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A21[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A22[i][n]);
                n+=1;
            }else {
                printf("\t  %f  ||",B2[i]);
            }
        }
    }
    //convierte en identidad A22
    copiaMatrizJ(AddInv, A22, dim-p, dim-p);
    inversa(AddInv, dim-p); //Se calcula la inversa de A dos x dos
    multiInverJ(A22, AddInv, dim-p, dim-p, dim-p, dim-p);
    //se realiza la sguiente operacion: B2 = A22^1 * B2
    productoMatVecJ(AddInv, B2, dim-p, dim-p);
    //se convierte en cero A12
    sumaProductJ(A12, A12copy, A22, p, dim-p, p, dim-p, dim-p, dim-p);
    //Se realiza la siguiente operacion: B1 = A12 * B2
    sumaProducVectorJ(B1, A12copy, B2, p, dim-p);
    //**************************    SIGUIENTE ITERACION ***************************************/
    printf("\n ITERACION 2 \n");
    for(i=0; i<p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A11[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A12[i][n]);
                n+=1;
            }else{
                printf("\t  %f  ||",B1[i]);
            }

        }
    }
    printf("\n");
    for(i=0; i<dim-p; i++){
        int n=0;
        printf("\n\n||");
        for(j=0; j<=dim; j++){
            if(j < p){
                printf("  %f  ",A21[i][j]);
            }else if(j>=p && j<dim){
                if(j == p){
                    printf("\t");
                }
                printf("  %f  ",A22[i][n]);
                n+=1;
            }else {
                printf("\t  %f  ||",B2[i]);
            }
        }
    }
    printf("\n\n Matriz A. ");
    mostrarMatriz(A, dim);
    printf("\n\n  Vector Solucion ");
    int h=0;
    for(i=0; i<dim; i++){
        if(i < p){
            Solucion[i] = B1[i];
        }else{
            Solucion[i] = B2[h];
            h += 1;
        }
    }
    printf("\n|-----------------------------|");
    for( i=0; i<dim; i++){
        printf("\n|   X%d   =    %.4f          |",i+1,Solucion[i]);
        printf("\n|-----------------------------|");
    }
    printf("\n\n Comprobacion ");
    printf("\n Matriz A x Vector solucion = Vector de terminos independintes. ");
    productoMatVecJ(A, Solucion, dim, dim);
    printf("\n A x Xs = B ");
    imprimeVector(Solucion, dim);

    //Liberacion de memoria
    liberarMatrizJ(AuuInv, p);
    liberarMatrizJ(AddInv, dim-p);
    liberarMatrizDdimJ(A12, p, dim-p);
    liberarMatrizJ(A11, p);
    liberaVectorJ(B1);
    liberarMatrizDdimJ(A21, dim-p, p);
    liberarMatrizDdimJ(A21copy, dim-p, p);
    liberarMatrizJ(A22, dim-p);
    liberaVectorJ(B2);
    liberarMatrizDdimJ(A12copy, p, dim-p);
    liberaVectorJ(Solucion);
}
int gauss(){
    MD A;
    pentero vector;
    int dim=0;
    printf("\n\t\tGAUSS - JORDAN PARTICIONADO ");
    printf("\n\n Ingresa el rango de la matriz: ");
    scanf("%d",&dim);
    A = generaMatrizJ(dim);
    leerMatrizJ(A, dim);
    printf("\n Matriz A ");
    mostrarMatrizJ(A, dim);
    int det = determinante(A, dim);
    if(det == 0){
        printf("\n El sistema no tiene solucion. ");
        return 13;
    }
    vector = generaVectorJ(dim);
    printf("\n Ingresa los elementos del vector de terminos independiente. ");
    ingresaVectorJ(vector, dim);
    printf("\n Vector de terminos independientes. ");
    imprimeVectorJ(vector, dim);
    gausJordanParticionadoJ(A, vector, dim);
    liberarMatrizJ(A, dim);
    liberaVectorJ(vector);
    printf("\n\n");
}

int portada(){
    printf("\nUniversidad Nacional Aut%cnoma de M%cxico.\n",162,130);
    printf("\nMatem%cticas Aplicadas y Computaci%cn.\n",160,162);
    printf("\nMetodos Numericos\n");
    printf("\nIntegrantes:\nMora Cano Paloma.\nGonzalez Pena Jose Osmar.\nJuarez Gallardo Isaac\n");
    system("pause");
    system("cls");
    return 0;
}

int main(){
	portada();
	int top;
	do{
	system("cls");
	int res,opc;
	printf("Escoja el metodo que desea realizar\n");
	printf("1.Metodo de biseccion\n");
	printf("2.Metodo de falsa posicion\n");
	printf("3.Metodo de Newton\n");
	printf("4.Metodo de la secante\n");
	printf("5.Inversion de matrices particionadas\n");
	printf("6.Gauss-Jordan particionado\n");
	printf("7.Metodo de intercambio\n");
	printf("8.Metodo de jacobi\n");
	printf("9.Metodo de Gauss-Seidel\n");
	printf("10.Metodo de Relajacion\n");
	printf("11.Metodo de Cholesky\n");
	printf("12.Metodo de Dolitlle\n");
	printf("13.Metodo de Crout\n");
	printf("opcion: ");
	scanf("%d",&opc);
	switch (opc){
		case 1:{
			do{
			
			system("cls");
			bis();
		printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
		scanf("%i",&res);
	}while(res==1);
	break;
}
		case 2:{
			do{
			system("cls");
			Fp();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
		}while(res==1);
			
		break;
	}
		case 3:{
			do{
			system("cls");
			newt();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
		case 4:{
			do{
			system("cls");
			sec();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
		case 5:{
			do{
			system("cls");
			Inver();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
		
		case 6:{
			do{
			system("cls");
			gauss();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
		
		case 7:{
			do{
			system("cls");
			Int();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
 		
 		case 8:{
 			do{
 		system("cls");
		jacobi();
		printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
		 break;
	}
		case 9:{
			{
			system("cls");
			seidel();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
		}
		
		case 10:{
			do{
			system("cls");
			rel();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
		break;
	}
		case 11:{
		do{
			system("cls");
			cholesky();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);		
			break;
		}
			
		case 12:{
			do{
			system("cls");
			dolittle();
			printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
			break;
	}
		case 13:{
			do{
		system("cls");
			Crou();
		printf("\n¿Desea repetir el metodo? (si=1 , no=0): ");
			scanf("%i",&res);
			}while(res==1);
		break;
}

}
printf("¿Desea regresar al menu principal? (si=1  no=0): ");
scanf("%i",&top);
}while(top==1);
system("pause");
return 0;
}

