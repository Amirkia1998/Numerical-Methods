#include <stdio.h>
#include <math.h>
#define COLEN 10         // size of coefficients array 
#define N 4             // number of rows of NxM matrix. can be changed from here
#define M 4             // number of columns of NxM matrix. can be changed from here
#define MAX_SIZE 20     // Maximum size of a matrix
#define nTable 5             // number of (x,y) in Newton-Gregory Interpolation 

// prototypes
void getCoefficients(double coefficients[], int *degree);
double calFunction(double coefficients[], int degree, double xInput);
void bisection(double coefficients[], int degree);
void regulaFalsi(double coefficients[], int degree);
void newRaph(double coefficients[], int degree);
double polDerivative(double coefficients[], int degree, double xInput);
void getCofactor(float A[N][N], float temp[N][N], int p, int q, int n);
void adjoint(float A[N][N],float adj[N][N]);
float determinant(float A[N][N], int n);
void inverse(float A[N][N], float inverse[N][N]);
void printMatrix(float A[N][N]);
void scanMatrix(float A[N][N]);
void gaussElimination(float upperTriangle[MAX_SIZE][MAX_SIZE], float gaussX[MAX_SIZE]);
void gaussSeidel(double equation[MAX_SIZE][MAX_SIZE], double x[MAX_SIZE]);
void numDif(double coefficients[], int degree);
double polIntegral(double coefficients[], int degree, double xInput);
void trapez(double coefficients[], int degree);
void simpson(double coefficients[], int degree);
float uCal(float u, int i);
int fact(i);
void gregoryNewton(float xTable[nTable], float yTable[nTable][nTable]);


//##########################################################################################################
//######################################## MAIN #########################################################
//##########################################################################################################
int main(void) {

    // variables
    int option;                                         // option is used for menu
    int i, j, degree;                                      // degree ---> degree of polynomia
    double coefficients[COLEN];                            // stores coefficients of a polynomia

    //-----------------------arrays for inverse of a matrix---------------------------------------------------- 
    float A[N][N];                              // Size of matrix (N) can be changed
    float adj[N][N];                            // To store adjoint of A[][]            
    float inv[N][N];                          // To store inverse of A[][]
    

    //------------------------------------arrays used for Gauss Elimination--------------------------------------------------------------------
    float upperTriangle[MAX_SIZE][MAX_SIZE];    // To store equations of Gauss Elimination method (it will be changed to upper triangle matrix)
    float gaussX[MAX_SIZE];                     // To store answers of linear equations solved by Gauss Elimination method


    //------------------------------------arays used for Gauss Seidel Elimination--------------------------------
    double equation[MAX_SIZE][MAX_SIZE];     // To store equations of Gauss Seidel 
    double x[MAX_SIZE];                     // To store initial values and next values of x, y, z
    

    //----------------------------Arrays used for Newton-Gregory Interpolation------------------------------------
    float xTable[nTable];        // x[] is used for x of (x,y)
    float yTable[nTable][nTable];     // y[][] is used for difference and y[i][0] used for y of (x,y)


    
    //--------------------------------------------------menu starts here------------------------------------------------
    do                  
    {       
        printf(" \n\n---------------------------------------------------------------------------------- \n");
        printf("1) Bisection Method \n2) Regula-Falsi Method \n3) Newton-Raphson \n4) Inverse of NxN matrix \n5) Gauss Elimination \n6) Gauss Seidel\n");
        printf("7) Numerical Differentation\n8) Trapezoidal method\n9) Simpson 1/3 method\n10) Gregory-Newton Interpolation\n0) Exit\n");
        printf("\nSelect one of the options above: ");
        scanf("%d", &option);
        printf(" \n-------------------------------------------------------------------------------------- \n");


        if (option == 1)    // bisection
        {
            getCoefficients(coefficients, &degree);
            bisection(coefficients, degree);
        }
        else if (option == 2)   // regula-falsi
        {
            getCoefficients(coefficients, &degree);
            regulaFalsi(coefficients, degree);
        }
        else if (option == 3)   // newton-raphson
        {
            getCoefficients(coefficients, &degree);
            newRaph(coefficients, degree);
        }
        else if (option == 4)   // matrix inverse
        {
            scanMatrix(A);
            printf("\nthe adjoint matrix is : \n");
            adjoint(A, adj);
            printMatrix(adj);
            printf("\ndeterminant is: ");
            printf("%f \n", determinant(A, N));

            inverse(A, inv);

            printf("\nInverse matrix is: \n");
            for ( i = 0; i < N; i++)
            {
                printf("| ");
                for ( j = 0; j < N; j++)
                {
                    printf("%f\t", inv[i][j]);
                }
                printf("\n");

            }
        }
        else if (option == 5)       // gauss elimination
        {
            gaussElimination(upperTriangle, gaussX);
        }
        else if (option == 6)       // gauss seidel
        {
            gaussSeidel(equation, x);
        }
        else if (option == 7)       // numerical differantiation
        {
            getCoefficients(coefficients, &degree);
            numDif(coefficients, degree);
        }
        else if (option == 8)       // trapezoidal method
        {
            getCoefficients(coefficients, &degree);
            trapez(coefficients, degree);
        }
        else if (option == 9)       // simpson method
        {
            getCoefficients(coefficients, &degree);
            simpson(coefficients, degree);
        }
        else if (option == 10)      // gregory newton method
        {
            gregoryNewton(xTable, yTable);
        }
        
        
        
        
    } while (option != 0);
    


}


//##########################################################################################################
//######################################## FUNCTIONS #########################################################
//##########################################################################################################

// function to get coefficients of a polynomia ---> ex) ax^2 + bx + c ---> a is 1. coefficients , b is 2. coefficients and ... 
void getCoefficients(double coefficients[], int *degree) {
    
    int i;

    do
    {
        printf("Please enter degree of polynomia (max 9): ");
        scanf("%d", degree);
    } while ( (*degree)<1 || (*degree)>=COLEN);

    for ( i = *degree; i >= 0; i--)
    {
        printf("enter %d. coefficient --->", (*degree)-i+1);
        scanf("%lf", &coefficients[i]);
    }
    
}


// calculate f(x) 
double calFunction(double coefficients[], int degree, double xInput) {

    double sum=0.0;
    int i;

    for ( i = 0; i < (degree+1); i++)
    {
        sum += coefficients[i] * pow(xInput,i);
    }

    return sum;
}


// bisection method
void bisection(double coefficients[], int degree){

    double epsilon;
    double a, b, c;
    
    do
    {
        printf("please enter epsilon: ");
        scanf("%lf", &epsilon);
    } while (epsilon <= 0);

    do
    {
        printf("enter a: ");
        scanf("%lf", &a);
        printf("enter b: ");
        scanf("%lf", &b);
    } while ( calFunction(coefficients, degree, a) * calFunction(coefficients, degree, b) >= 0 );

    printf("\n");
 
    c = a;
 
    while ((b-a) >= epsilon)
    {
        c = (a+b)/2;
        //if (calFunction(coefficients, degree ,c) == 0.0){
        //    printf("Root = %lf\n",c);
        //    break;
        //}
        if (calFunction(coefficients, degree ,c) * calFunction(coefficients, degree ,a) < 0){
                printf("Root = %lf\n",c);
                b = c;
        }
        else{
                printf("Root = %lf\n",c);
                a = c;
        }
    }
}


// regula-falsi method
void regulaFalsi(double coefficients[], int degree) {

    double a, b, c;
    double epsilon, error=10;
    int iteration=1;

    do
    {
        printf("please enter epsilon: ");
        scanf("%lf", &epsilon);
    } while (epsilon <= 0);

    //do
    //{
        printf("enter a: ");
        scanf("%lf", &a);
        printf("enter b: ");
        scanf("%lf", &b);
        printf("\n");
    //} while ( calFunction(coefficients, degree, a) * calFunction(coefficients, degree, b) >= 0 );

    while (error > epsilon)
    {
        
        c = ( a*calFunction(coefficients, degree, b) - b*calFunction(coefficients, degree, a) )
        /( calFunction(coefficients, degree, b) - calFunction(coefficients, degree, a) );
        
        iteration++;
        error = ( b-c )/( 1.0*pow(2,iteration) );

        if (calFunction(coefficients, degree, a) * calFunction(coefficients, degree, c) < 0)
        {
            printf("Root = %lf\n", c);
            b = c;
        }
        else
        {
            printf("Root = %lf\n", c);
            a = c;
        }
    }
}


// newton-raphson ,ethod
void newRaph(double coefficients[], int degree) {
    double x0, x1;
    double epsilon, error;
    double guess;
    int iter=0;

    printf("enter a: ");
    scanf("%lf", &guess);
    printf("enter epsilon: ");
    scanf("%lf", &epsilon);
    printf("\n");

    x0 = guess;
    do
    {
        iter++;
        x1 = x0 - (calFunction(coefficients, degree, x0) / polDerivative(coefficients, degree, x0));
        error = abs(x1-x0);
        x0 = x1;
        printf("Iteration:%d   Root=%lf    Error=%lf \n",iter, x1, error);
    } while (error > epsilon && iter!=100);
    

    if (iter == 100)        // maximum number of iterations is limited to 100
    {
        printf("\nIt might be convergant !!! Newton-Raphson FAILED :) ");
    }

}


// calculates derivate of polynomia
double polDerivative(double coefficients[], int degree, double xInput) {
    double sum=0.0;
    int i;
    int coDer[COLEN];       // coefficients of df(x)

    for ( i = degree; i > 0; i--)
    {
        coDer[i-1] = coefficients[i]*i;
    }

    for ( i = 0; i < (degree); i++)
    {
        sum += coDer[i] * pow(xInput,i);
    }

    return sum;
}


//----------------------- functions used for calculating inverse of NxN matrix START here-------------------------
void getCofactor(float A[N][N], float temp[N][N], int p, int q, int n)
{
    int row, col;
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (row = 0; row < n; row++)
    {
        for (col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];
  
                // Row is filled, so increase row index and reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
  

// Recursive function for finding determinant of matrix. n is current dimension of A[][].
float determinant(float A[N][N], int n)
{
    float D = 0.0; // Initialize result
    float sign = 1.0;  // To store sign multiplier
    float temp[N][N]; // To store cofactors
    int f;
  
    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];
  
     // Iterate for each element of first row
    for (f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);
  
        // terms are to be added with alternate sign
        sign = -sign;
    }
  
    return D;
}
  

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(float A[N][N],float adj[N][N])
{
    // temp is used to store cofactors of A[][]
    int sign = 1;
    float temp[N][N];
    int i, j;

    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }
  
    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);
  
            // sign of adj[j][i] positive if sum of row and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
  
            // Interchanging rows and columns to get the transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}


// Function to calculate and store inverse
void inverse(float A[N][N], float inverse[N][N])
{
    float det = determinant(A, N);
    float adj[N][N];
    int i, j;
    
    // Find determinant of A[][]
    
  
    // Find adjoint
    adjoint(A, adj);
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/(1.0*det);

}


//----------------------- functions used for calculating inverse of NxN matrix END here-------------------------

void printMatrix(float A[N][N]) {
    int i, j;

    for ( i = 0; i < N; i++)
    {   
        printf("| ");
        for ( j = 0; j < N; j++)
        {
            printf("%f\t", A[i][j]);
        }
        
        printf("\n");
        
    }
    
}

void scanMatrix(float A[N][N]) {
    int i, j;

    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            printf("enter %dx%d element: ", i+1, j+1);
            scanf("%f", &A[i][j]);
        }
        
    }
    
}


// Gauss elimination 
void gaussElimination(float upperTriangle[MAX_SIZE][MAX_SIZE], float gaussX[MAX_SIZE]) {
    int i, j, k;
    float co, sum=0.0;

    for ( i = 1; i <= N; i++)
    {
        for ( j = 1; j <= N+1; j++)
        {
            printf("enter %dx%d element: ", i, j);
            scanf("%f", &upperTriangle[i][j]);
        }
        
    }

    for(j=1; j<=N; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=1; i<=N; i++)
        {
            if(i>j)
            {
                co = upperTriangle[i][j] / upperTriangle[j][j];
                for(k=1; k<=N+1; k++)
                {
                    upperTriangle[i][k] = upperTriangle[i][k ]- co*upperTriangle[j][k];
                }
            }
        }
    }

    printf("\nUpper-Triangular matrix:\n");
    for ( i = 1; i <= N; i++)
    {   
        printf("|");
        for ( j = 1; j <= N+1; j++)
        {
            printf("%f\t", upperTriangle[i][j]);
            if (j == N)                         // for aesthetic
            {
                printf("|\t");
            }
            
        }
        printf("\n");
    }
    
    gaussX[N]=upperTriangle[N][N+1]/upperTriangle[N][N];

    for(i=N-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=N; j++)
        {
            sum += (upperTriangle[i][j]) * (gaussX[j]);
        }
        gaussX[i]=(upperTriangle[i][N+1]-sum)/upperTriangle[i][i];
    }
    printf("\nThe solution is: \n");
    for(i=1; i<=N; i++)
    {
        printf("\nx%d=%f\t",i,gaussX[i]); /* x1, x2, x3 are the required solutions*/
    }
    
}
    

// gauss seidel    
void gaussSeidel(double equation[MAX_SIZE][MAX_SIZE], double x[MAX_SIZE]) {
    int n, i, j, k, flag=0, count=0;
    double epsilon, y, tmp;

    printf("enter number of equations: ");        
    scanf("%d", &n);                    //Input number of equations

    printf("enter elements of matrix:\n");
    for (i=0;i<n;i++) 
    {
        for (j=0;j<n+1;j++) 
        {
            printf("%dx%d: ", i+1, j+1);
            scanf("%lf", &equation[i][j]);      // Input co-efficients of equations
        }
    }

    printf("\nenter initial values of variables:\n");   // initialize x, y, z
    for (i=0;i<n;i++) 
    {
        printf("enter x%d: ", i+1);
        scanf("%lf", &x[i]);            // x is for values of x, y, z (initial and next values)
    }

    printf("\nEnter the epsilon:\n");
    scanf("%lf", &epsilon);             // input epsilon

    for (i=0;i<n;i++)   //Pivotisation to make the equations diagonally dominant    
    {                
        for (k=i+1;k<n;k++)
        {
            if (abs(equation[i][i])<abs(equation[k][i]))
            {
                for (j=0;j<=n;j++)
                {
                    tmp=equation[i][j];
                    equation[i][j]=equation[k][j];
                    equation[k][j]=tmp;
                }
            }
        }
    }

    printf("\nIter \t");      // for aesthetic
    for(i=0;i<n;i++)            // for aesthetic
        printf("x%d \t\t\t", i+1);      // for aesthetic
    printf("\n----------------------------------------------------------------------------------"); // for aesthetic

    do                            //Perform iterations to calculate x1,x2,...xn
    {
        printf("\n\t");
        for (i=0;i<n;i++)                //Loop that calculates x1,x2,...xn
        {
            y=x[i];
            x[i]=equation[i][n];

            for (j=0;j<n;j++)
            {
                if (j!=i)
                x[i]=x[i]-equation[i][j]*x[j];
            }

            x[i]=x[i]/equation[i][i];
            if (abs(x[i]-y)<=epsilon)            //Compare the new value with the last value
                flag++;
            printf("%lf \t\t", x[i]);
        }
        printf("\n");
        count++;   
    }while(flag<n);     //If the values of all the variables don't differ from their previious values 
                        //with error more than epsilon then flag must be n and hence stop the loop

    printf("\nFinal Solution:  ");
    for ( i = 0; i < n; i++)
    {
        printf("x%d=%lf \t", i+1, x[i]);
    }
    
}

    
// numerical differentiation
void numDif(double coefficients[], int degree){
    
    double h, x, df, exact;
    int option;

    printf("\nplease enter x: ");
    scanf("%lf", &x);
    printf("please enter h: ");
    scanf("%lf", &h);

    exact = polDerivative(coefficients, degree, x);
    
    do
    {
        printf("\n1) Forward differentiation\n2) Backward differentiation\n3) Central differentiation\n");
        printf("Choose one of the options above: ");
        scanf("%d", &option);

        if (option == 1)    // forward
        {
            df = ( calFunction(coefficients, degree, x+h) - calFunction(coefficients, degree, x) ) / h;
            printf("\ndf = %lf\n", df);
        }
        else if (option == 2)       // backward
        {
            df = ( calFunction(coefficients, degree, x) - calFunction(coefficients, degree, x-h) ) / h;   
            printf("\ndf = %lf\n", df);
        }
        else if (option == 3)       // central
        {
            df = ( calFunction(coefficients, degree, x+h) - calFunction(coefficients, degree, x-h) ) / (2*h);
            printf("\ndf = %lf\n", df);
        }

    } while (option!=1 && option!=2 && option!=3);

    printf("Real df = %lf\n", exact);
    printf("Error = %lf", df-exact);
    
}


// calculates exact integral of a polinom
double polIntegral(double coefficients[], int degree, double xInput) {
    
    double coInt[COLEN+1] = {0};        // note: always coInt[0] = 0 beacause we dont know exact value of C after getting integral so we assume it as 0
    double sum = 0.0;
    int i;
    
    for ( i = degree; i >= 0; i--)
    {
        coInt[i+1] = (coefficients[i] / (i+1)); 
    }

    for ( i = 0; i < (degree+2); i++)         // i could be strated from 1 as we know coInt[0] is always 0;
    {
        sum = sum + coInt[i] * pow(xInput,i);
    }

    return sum;
}


// trapez method
void trapez(double coefficients[], int degree) {

    int n, k=1, upper, lower;         // lower---> X0     upper---> Xn
    double h, s1, s2=0, st;
    double exact, newX;

    printf("enter n: ");
    scanf("%d", &n);
    printf("enter lower interval: ");
    scanf("%d", &lower);
    printf("enter upper interval: ");
    scanf("%d", &upper);

    exact =  polIntegral(coefficients, degree, 1.0*upper) - polIntegral(coefficients, degree, 1.0*lower); 
    h = (1.0*upper-1.0*lower)/n;
    printf("\nh = %lf\n\n", h);

    s1 = h * ( calFunction(coefficients, degree, upper) + calFunction(coefficients, degree, lower) ) / 2;
        
    while (k<n)
    {
        newX = lower + (1.0*k)*h;
        printf("x%d = %lf \t", k ,newX);
        s2 += calFunction(coefficients, degree, newX);
        printf("y%d = %lf\n", k, s2);
        k++;
    }
        
    st = s1 + h*s2;
    
    printf("\nCalculated Integral = %lf \t Exact Integral = %lf \t Error = %lf", st, exact, st-exact);
}


// simpson method
void simpson(double coefficients[], int degree) {

    double upper, lower, h, s1=0, s2=0, s3=0, st=0, newX, exact;   // lower---> X0     upper---> Xn
    int n, k;

    printf("enter n: ");
    scanf("%d", &n);
    printf("enter lower interval: ");
    scanf("%lf", &lower);
    printf("enter upper interval: ");
    scanf("%lf", &upper);

    exact =  polIntegral(coefficients, degree, upper) - polIntegral(coefficients, degree, lower);
    h = (upper-lower)/n;
    printf("\nh = %lf\n\n", h);

    s1 = ( calFunction(coefficients, degree, upper) + calFunction(coefficients, degree, lower) );

    for ( k = 1; k < n; k=k+2)      // k = 1,3,5,...
    {
        newX = lower + k*h;
        printf("x%d = %lf \t", k ,newX);
        s2 += calFunction(coefficients, degree, newX);
        printf("y%d = %lf\n", k, s2);
    }
    printf("\n");
    for ( k = 2; k < n-1; k=k+2)        // k = 2,4,6,...
    {
        newX = lower + k*h;
        printf("x%d = %lf \t", k ,newX);
        s3 += calFunction(coefficients, degree, newX);
        printf("y%d = %lf\n", k, s3);
    }
    
    st = (h/3.0)*(s1 + 4.0*s2 + 2.0*s3);

    printf("\nCalculated Integral = %lf \t Exact Integral = %lf \t Error = %lf", st, exact, st-exact);
}


// in gregory interpolation, u is ---> x0(x0-x1)(x0-x2)....
float uCal(float u, int i)
{
    int j;
    float temp = u;
    for (j = 1; j < i; j++)
        temp = temp * (u + j);
    return temp;
}


// calculates factoriel of input i
int fact(i)
{
    int j;
    int f = 1;
    for ( j = 2; j <= i; j++)
        f *= j;
    return f;
}


void gregoryNewton(float xTable[nTable], float yTable[nTable][nTable]) {

    int i, j;
    float value;    // value of x in F(x)
    float sum, u;

    for ( i = 0; i < nTable; i++)        // input x
    {
        printf("please enter x%d: ", i);
        scanf("%f", &xTable[i]);
    }

    for ( i = 0; i < nTable; i++)        // input y
    {
        printf("please enter y%d: ", i);
        scanf("%f", &yTable[i][0]);
    }

    printf("please enter value of x in F(x): ");
    scanf("\n%f", &value);

    // Calculating the backward difference table
    for (i = 1; i < nTable; i++) {
        for (j = nTable - 1; j >= i; j--)
            yTable[j][i] = yTable[j][i - 1] - yTable[j - 1][i - 1];
    }

    printf("\n");
    // Displaying the backward difference table
    for ( i = 0; i < nTable; i++) {
        for ( j = 0; j <= i; j++)
            printf("%f \t", yTable[i][j]);

        printf("\n");                   
    }

    // Initializing u and sum
    sum = yTable[nTable - 1][0];
    u = (value - xTable[nTable - 1]) / (xTable[1] - xTable[0]);
    for ( i = 1; i < nTable; i++) {
        sum = sum + (uCal(u, i) * yTable[nTable - 1][i]) / (1.0*fact(i));
    }
 
    printf("\n\nF(%f) = %f \n", value, sum);
    
}

