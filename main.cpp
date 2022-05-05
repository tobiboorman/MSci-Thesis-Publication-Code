//  Quantum Circuit Program
//
//  Authors:
//
//      -Tobias Boorman
//      -Dr Marcin Szyniszewski
//
//  Inputs for shell script:
//
//      -circuitLength   (Number of qubits)
//      -piece           (See: How mutual information (MI) works)
//      -measProbability (Probability for measurement, program scales this with dt)
//      -measStrength    (Measurement Strength)
//      -xiStatic        (Strength of static disorder)
//      -xiTemporal      (Strength of uniform temporal disorder, program scales this with dt)
//      -xiNonstatic     (Strength of non-static random disorder, program scales this with dt)
//      -dt              (Time step increment)
//      -startReadings   (When to start entropy readings - real time)
//      -entropyTiming   (Timing between entropy readings in units of dt) older versions used units of real time but this led to problematic floating point errors
//      -simTime         (Time over which sim is performed - real time)

// Mutual information calculation
//
//  MI is calculated assuming the system is split like this:
//
//      A            C            B            D
//  |=======|=================|=======|=================|
//   L/piece                   L/piece
//
//  Then MI = S(A) + S(B) - S(A u B).
//
//  If piece=4, then tripartite mutual information (TMI) is calculated as
//  TMI = S(A) + S(B) + S(C) - S(A u B) - S(A u C) - S(B u C) + S(A u B u C)

// ----- Preamble ----- //

//Packages:
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"
#include "Eigen/Dense"
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <cmath>
#include <random>
#include <chrono>

//Name spaces:
using namespace Eigen;
using namespace std;

// ----- Objects ----- //

//Random number generator:
mt19937_64 generator(chrono::system_clock::now().time_since_epoch().count());

//Misc definition of real and imaginary unity (double):
complex<double> Id(0.0, 1.0), Ud(1.0,0.0);

// ----- Distributions ----- //

//Uniform distribution used in measurement process:
uniform_real_distribution<double> random_real(0.0,1.0);

// ------ Functions ----- //

//Integer exponentiation (Taken directly from https://stackoverflow.com/a/101613, speeds up the integer exponentiation process):
int ipow(int base, int exp)
{
    int result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

//The eigensolver for exponentiation of 4x4 matrices (with i and dt):
SelfAdjointEigenSolver< Matrix<complex<double>,4,4> > eigensolver4(4);

//Matrix exponent exp(-i dt H), where H must be Hermitian:
Matrix<complex<double>,4,4> matrixExp(Matrix<complex<double>,4,4> H, double t)
{
    eigensolver4.compute(H);
    return eigensolver4.eigenvectors() * (-complex<double>(0.0,t)*eigensolver4.eigenvalues().array()).exp().matrix().asDiagonal() * eigensolver4.eigenvectors().adjoint();
}

//Kronecker product of two matrices or vectors:
MatrixXcd Kronecker(MatrixXcd A,MatrixXcd B){

    int nrowa = A.rows();
    int nrowb = B.rows();
    int ncola = A.cols();
    int ncolb = B.cols();
    div_t k_loc;

    //C will be the output of this function:
    MatrixXcd C(nrowa*nrowb,ncola*ncolb);

    //Process for performing a Kronecker product:

    //Over the elements of A:
    for(int k=0;k < A.size();k++){

        //The quotient and remainder, used for identifying position of element in A in terms of k:
        k_loc = div(k,ncola);

        //Over the rows of B:
        for(int i=0;i<nrowb;i++){
            //Over the columns of B:
            for(int j=0;j<ncolb;j++){
                //Define C element wise:
                C(k_loc.quot*nrowa + i, k_loc.rem*ncola + j) = A(k_loc.quot*(nrowa-1),k_loc.rem*(ncola-1))*B(i,j);
            }
        }
    }
    return C;
}

// ----- Classes ----- //

//The quantum circuit class:
class circuit{

    //Access Specifier:
    public:

    //The constructor:
    circuit (int, int);

    //Characteristic lengths ...

    //... Relating to bipartite entanglement entropy:
    int length; //Length of circuit
    int stateLength; //The length of the state vector
    int halfstateLength; //The length of the state for half the number of qubits

    //... Relating to MI:
    int piece;
    int lenA;
    int lenApr;
    int lenSub;
    int lenSubpr;
    int lenAA;

    //The circuit state vector:
    complex<double>* state;

    //Density matrices for entropy and MI calculations:
    MatrixXcd rhoBpe;
    MatrixXcd rhoAIaux;
    MatrixXcd rhoAIauxlong;
    MatrixXcd rhoA;
    MatrixXcd rhoB;
    MatrixXcd rhoAB;

    //Eigensolvers for entropy and MI readings:
    SelfAdjointEigenSolver<MatrixXcd> eigensolverBptent;
    SelfAdjointEigenSolver<MatrixXcd> eigensolverA;
    SelfAdjointEigenSolver<MatrixXcd> eigensolverAB;

    //Variables for entropy and MI calculations:
    double mutualinf;
    double S_A;
    double S_B;
    double S_C;
    double S_AB;
    double S_AC;
    double S_BC;
    double S_ABC;
    double TMI;

    //Set the system to the Neel State:
    void neelState();

    //Print the state of the system:
    void giveState(){
        cout <<  state << endl;
    }

    //Give the state norm:
    void giveNorm(){
        double sum = 0;
        for(int i = 0; i < stateLength; i++)
            sum += norm(state[i]);
        cout << sqrt(sum) << endl;

    }

    //Normalise the state:
    void normState(){
        double sum = 0;
        for(int i = 0; i < stateLength; i++)
            sum += norm(state[i]);
        sum = sqrt(sum);
        for(int i = 0; i < stateLength; i++)
            state[i] /= sum;
    }

    //Apply double-site process to sites j and j+1:
    void applySites(int,Matrix<complex<double>,4,4>&);

    //Parameters and distributions for the measurement process:
    double lambda;
    double delta;
    normal_distribution<double> Gauss_U, Gauss_D;

    //Set up the measurement process for the circuit:
    void measurementParameters(double, double);

    //Perform a generalised measurement of the jth qubit:
    void measureState(int);

    //Determine the entropy of system:
    double findEntropy(SelfAdjointEigenSolver< MatrixXcd > &, MatrixXcd &, int);

    //Determine the MI:
    void findMI();
};

//Circuit constructor:
circuit::circuit (int x, int p){

    //System length:
    length = x;

    //Piece variable for MI:
    piece = p;

    //Set characteristic scales:
    lenA = ipow(2,x/p);
    lenApr = ipow(2,x-(x/p));
    lenSub = ipow(2,x/2-(x/p));
    lenSubpr = ipow(2,x/2+(x/p));
    lenAA = lenA*lenA;
    stateLength = ipow(2,x);
    halfstateLength = ipow(2,x/2);

    //Initialise system state and density matrices:
    state = new complex<double>[stateLength];
    rhoBpe = MatrixXcd::Zero(halfstateLength,halfstateLength);
    rhoAIaux = MatrixXcd::Zero(lenApr,lenA);
    rhoAIauxlong = MatrixXcd::Zero(lenSub*lenSub,lenA*lenA);
    rhoA = MatrixXcd::Zero(lenA,lenA);
    rhoB = MatrixXcd::Zero(lenA,lenA);
    rhoAB = MatrixXcd::Zero(lenA*lenA,lenA*lenA);

    //Set up eigensolvers for entropy and MI calculations:
    SelfAdjointEigenSolver<MatrixXcd> eigensolverBptent(halfstateLength);
    SelfAdjointEigenSolver<MatrixXcd> eigensolverA(lenA);
    SelfAdjointEigenSolver<MatrixXcd> eigensolverAB(lenA*lenA);
}

//Set the system to the anti-ferromagnetic Neel State:
void circuit::neelState (){

    int x = length;

    //index of unity in Neel State:
    int index = 0;

    //Find the index:
    while(x>=2){index += ipow(2,x-2); x -= 2;}
    state = new complex<double>[stateLength];
    state[index] = 1.0;
}

void circuit::applySites(int j, Matrix<complex<double>,4,4>& A){

    //Temporary variables:
    complex <double> t1,t2,t3,t4;

    //State length:
    if(j < length - 1)
    {
        //Number of blocks:
        int NB = ipow(2,j);

        //Block width:
        int BW = stateLength/NB;

        //Quarter block width:
        int QBW = BW/4;

        //Block starting location:
        int BSL;

        //Subblock starting location:
        int SSL;

        //Loop over blocks:
        for(int i = 0; i < NB; i ++){
            BSL = i*BW;
            for(int k = 0; k < QBW; k++){

                SSL = BSL+k;

                //Temporary variables:
                t1 = state[SSL      ];
                t2 = state[SSL+  QBW];
                t3 = state[SSL+2*QBW];
                t4 = state[SSL+3*QBW];

                //Process itself:
                state[SSL      ] = A(0,0)*t1 + A(0,1)*t2 + A(0,2)*t3 + A(0,3)*t4;
                state[SSL+  QBW] = A(1,0)*t1 + A(1,1)*t2 + A(1,2)*t3 + A(1,3)*t4;
                state[SSL+2*QBW] = A(2,0)*t1 + A(2,1)*t2 + A(2,2)*t3 + A(2,3)*t4;
                state[SSL+3*QBW] = A(3,0)*t1 + A(3,1)*t2 + A(3,2)*t3 + A(3,3)*t4;
            }
        }
    }
    else
    {
        //Last spin
        int LS = stateLength/2;

        for(int p = 0; p < LS; p+=2){

            //Temporary variables:
            t1 = state[p];
            t2 = state[LS+p];
            t3 = state[p+1];
            t4 = state[LS+p+1];

            //Process itself:
            state[p]      = A(0,0)*t1 + A(0,1)*t2 + A(0,2)*t3 + A(0,3)*t4;
            state[LS+p]   = A(1,0)*t1 + A(1,1)*t2 + A(1,2)*t3 + A(1,3)*t4;
            state[p+1]    = A(2,0)*t1 + A(2,1)*t2 + A(2,2)*t3 + A(2,3)*t4;
            state[LS+p+1] = A(3,0)*t1 + A(3,1)*t2 + A(3,2)*t3 + A(3,3)*t4;
        }
    }
}

//Set up measurement parameters for an instance of the circuit class:
void circuit::measurementParameters(double l, double d){

    //Measurement strength and width of distribution:
    lambda = l;
    delta = d;

    //Set the distributions:
    Gauss_U = normal_distribution<double>(-lambda,delta/sqrt(2.0));
    Gauss_D = normal_distribution<double>(lambda,delta/sqrt(2.0));
}

//Perform a generalised measurement of jth qubit:
void circuit::measureState(int j){

    //Step 1 - Variables that store the sum of the norms of the spin up and spin down components respectively:
    double S_U, S_D;
    S_U = S_D = 0;

    //Step 2 - Start adding summands to S_U and S_D:
    int sizeB = ipow(2,j), sizeA = stateLength/sizeB/2, sizeB2 = sizeB*2, temp, temp2;
    for(int a = 0; a < sizeA; a++)
    {
        temp = a*sizeB2;
        for(int b = 0; b < sizeB; b++)
            S_U += norm(state[temp + b]);
    }
    S_D = 1.0 - S_U; //S_U is the probability of having spin up

    //Step 3 - Sample x:
    double x = (random_real(generator) < S_U ? Gauss_U(generator) : Gauss_D(generator));

    //Step 5 - Update the state:
    double cU = 1.0/sqrt( S_U + S_D * exp(+4.0*x*lambda) );
    double cD = 1.0/sqrt( S_D + S_U * exp(-4.0*x*lambda) );
    for(int a = 0; a < sizeA; a++)
    {
        temp = a*sizeB2;
        temp2= a*sizeB2 + sizeB;
        for(int b = 0; b < sizeB; b++)
            state[temp  + b] *= cU;
        for(int b = 0; b < sizeB; b++)
            state[temp2 + b] *= cD;
    }

}

//Diagonalise the reduced density matrix and return entropies:
double circuit::findEntropy(SelfAdjointEigenSolver< MatrixXcd > & eigensolver, MatrixXcd & rho, int size)
{
    eigensolver.compute(rho,EigenvaluesOnly);
    double entropy  = 0, eig;
    for(int i = 0; i < size; i++)
    {
        eig = real(eigensolver.eigenvalues()(i));
        if(eig>1.e-15)
        {
            entropy  -= eig * log(eig);
        }
    }
    return entropy;
}

//Calculate the MI (output: bipartite entropy, MI, tripartite mutual information):
void circuit::findMI()
{
    //Set variables to zero:
    mutualinf = 0, S_A = 0, S_B = 0, S_C = 0, S_AB = 0, S_AC = 0, S_BC = 0, S_ABC = 0, TMI = 0;

    //Find the bipartite entropy, S_AC (halfstateLength, halfstateLength):
    for(int b1 = 0; b1 < halfstateLength; b1++)
    for(int b2 = 0; b2 < halfstateLength; b2++)
        rhoBpe(b2,b1) = state[b1*halfstateLength + b2];
    rhoBpe = rhoBpe * rhoBpe.adjoint();

    S_AC = findEntropy(eigensolverBptent,rhoBpe, halfstateLength);

    cout << S_AC;

    if(length%piece==0)
    {

        // S_A (lenA, lenApr)
        for(int b1 = 0; b1 < lenA; b1++)
        for(int b5 = 0; b5 < lenApr; b5++)
            rhoAIaux(b5,b1) = state[b1*lenApr + b5];
        rhoA = rhoAIaux.adjoint()*rhoAIaux;

        S_A = findEntropy(eigensolverA, rhoA, lenA);

        // S_(AuB) (lenA, lenSub, lenA, lenSub)
        for(int b1 = 0; b1 < lenA; b1++)
        for(int b3 = 0; b3 < lenA; b3++)
        for(int b5 = 0; b5 < lenSub; b5++)
        for(int b6 = 0; b6 < lenSub; b6++)
            rhoAIauxlong(b5*lenSub+b6,b1*lenA+b3) = state[b1*lenApr + b3*lenSub + b5*halfstateLength + b6];
        rhoAB = rhoAIauxlong.adjoint()*rhoAIauxlong;

        S_AB = findEntropy(eigensolverAB, rhoAB, lenAA);

        // S_B (halfstateLength, lenA, lenSub)
        for(int b1 = 0; b1 < lenA; b1++)
        for(int b5 = 0; b5 < halfstateLength; b5++)
        for(int b6 = 0; b6 < lenSub; b6++)
            rhoAIaux(b5*lenSub + b6, b1) = state[b5*halfstateLength + b1*lenSub + b6];
        rhoA = rhoAIaux.adjoint()*rhoAIaux;

        S_B = findEntropy(eigensolverA, rhoA, lenA);

        mutualinf = S_A - S_AB + S_B;

        cout << "," << mutualinf;

        if(piece==4) // Calculating the rest of entropies for tripart MI
        {
            // S_C (lenA, lenSub, halfstateLength)
            for(int b1 = 0; b1 < lenSub; b1++)
            for(int b5 = 0; b5 < lenA; b5++)
            for(int b6 = 0; b6 < halfstateLength; b6++)
                rhoAIaux(b5*halfstateLength + b6, b1) = state[b5*lenApr + b1*halfstateLength + b6];
            rhoA = rhoAIaux.adjoint()*rhoAIaux;

            S_C = findEntropy(eigensolverA, rhoA, lenSub);

            // S_BC (lenA, halfstateLength, lenSub)
            for(int b1 = 0; b1 < halfstateLength; b1++)
            for(int b5 = 0; b5 < lenA; b5++)
            for(int b6 = 0; b6 < lenSub; b6++)
                rhoAIauxlong(b5*lenSub + b6, b1) = state[b5*lenApr + b1*lenSub + b6];
            rhoAB = rhoAIauxlong.adjoint()*rhoAIauxlong;

            S_BC = findEntropy(eigensolverAB, rhoAB, halfstateLength);

            // S_ABC = S_D (lenSubpr, lenSub)
            for(int b1 = 0; b1 < lenSub; b1++)
            for(int b5 = 0; b5 < lenSubpr; b5++)
                rhoAIaux(b5,b1) = state[b5*lenSub + b1];
            rhoA = rhoAIaux.adjoint()*rhoAIaux;

            S_ABC = findEntropy(eigensolverA, rhoA, lenSub);

            TMI = S_A + S_B + S_C - S_AB - S_AC - S_BC + S_ABC;

            cout << "," << TMI;
        }
    }
    cout << endl;
}

// ----- Main Body ----- //

int main(int argc, char* argv[]){

    //Set a large output precision:
    cout << fixed;
    cout.precision(15);

    //Default variables:
    int circuitLength = 10; //Number of qubits
    double Jx, Jy, Jz; //Spin couplings
    Jx = 1;
    Jy = 1;
    Jz = 1;
    double dt = 0.01; //The real time per step
    int numSteps = 100; //The number of steps in total
    double measProbability; //Measurement frequency (0-1)
    double measStrength, measStrength0 = 0; //Measurement strength (user input, scaled user input)
    double stdDev = 1; //Standard deviation (measurement process)
    double xiStatic = 0; //Static disorder strength
    double xiNonstatic, xiNonstatic0 = 0; //Non-static disorder strength (user input, scaled Hamiltonian input)
    double xiTemporal, xiTemporal0 = 0; //Temporal disorder strength (user input, scaled Hamiltonian input)
    int entropyTiming = 1; //Steps between entropy readings (In units of dt)
    double startReadings = 1; //When to start the entropy readings (In real time)
    double simTime; //The real time (numSteps*dt) over which the simulation is performed
    double piece = 4; //piece for MI

    //Shell inputs:
    if(argc == 12){

        //Inputs:
        circuitLength = atoi(argv[1]);
        piece = atoi(argv[2]);
        measProbability = atof(argv[3]);
        measStrength = atof(argv[4]);
        xiStatic = atof(argv[5]);
        xiTemporal = atof(argv[6]);
        xiNonstatic = atof(argv[7]);
        dt = atof(argv[8]);
        startReadings = atof(argv[9]);
        entropyTiming = atoi(argv[10]);
        simTime = atof(argv[11]);

        //Set the number of steps:
        numSteps = simTime/dt;

        //Scale the disorder strengths, for the Hamiltonian input:
        xiTemporal0 = xiTemporal/sqrt(dt);
        xiNonstatic0 = xiNonstatic/sqrt(dt);

        //Scale the measurement strength, for measureState input (caution: this scaling has only been confirmed for measProbability = 1, further work needed to understand interplay):
        measStrength0 = measStrength*sqrt(dt);

    }

    //Distributions:
    uniform_real_distribution<double> disorderDistribution(-1,1);

    //The Pauli matrices and their Kronecker squared value, useful for the Suzuki-Trotter decomposition:
    Matrix2cd p_0;
    Matrix2cd p_x;
    Matrix2cd p_y;
    Matrix2cd p_z;

    p_0 << Ud,0.0*Id,
           0.0*Id,Ud;

    p_x << 0.0*Id,1.0*Ud,
           1.0*Ud,0.0*Id;

    p_y << 0.0*Id, -1.0*Id,
           1.0*Id,0.0*Id;

    p_z << 1.0*Ud,0.0*Id,
           0.0*Id,-1.0*Ud;

    MatrixX4cd p_0_2 = Kronecker(p_0,p_0);
    MatrixX4cd p_x_2 = Kronecker(p_x,p_x);
    MatrixX4cd p_y_2 = Kronecker(p_y,p_y);
    MatrixX4cd p_z_2 = Kronecker(p_z,p_z);

    //Matrix for two-site interaction layer:
    Matrix4cd couplingTerm = Jx*p_x_2+Jy*p_y_2+Jz*p_z_2;

    //Static disorder vectors:
    VectorXd h_s_x = VectorXd::Zero(circuitLength);
    VectorXd h_s_y = VectorXd::Zero(circuitLength);
    VectorXd h_s_z = VectorXd::Zero(circuitLength);

    //Fill out static disorder vectors:
    for(int i = 0; i < circuitLength; i++){
        h_s_x[i] = disorderDistribution(generator);
        h_s_y[i] = disorderDistribution(generator);
        h_s_z[i] = disorderDistribution(generator);
    }

    //Store and fill out the static disorder Hamiltonians:
    MatrixX2cd staticH[circuitLength];
    for(int j = 0; j < circuitLength; j++){
        staticH[j] = xiStatic*(h_s_x[j]*p_x + h_s_y[j]*p_y + h_s_z[j]*p_z);
    }

    //Store the Suzuki-Trotter matrices:
    Matrix<complex<double>,4,4> stMatrices[circuitLength];

    //Store the non-static Hamiltonians:
    MatrixX2cd nonstaticH[circuitLength];

    //Initialise state and measurement process:
    circuit QC1(circuitLength, piece);
    QC1.neelState();
    QC1.measurementParameters(measStrength0, stdDev);

    double h_t_x, h_t_y, h_t_z;
    double h_r_x, h_r_y, h_r_z;

    //Experiment:
    for(int i = 0; i < numSteps; i++){

        //Temporal disorder Hamiltonian:
        h_t_x = disorderDistribution(generator);
        h_t_y = disorderDistribution(generator);
        h_t_z = disorderDistribution(generator);
        MatrixX2cd temporalH = xiTemporal0*(h_t_x*p_x + h_t_y*p_y + h_t_z*p_z);

        //Fill out non-static term:
        for(int j = 0; j < circuitLength; j++){

            //Fill out non-static disorder terms:
            h_r_x = disorderDistribution(generator);
            h_r_y = disorderDistribution(generator);
            h_r_z = disorderDistribution(generator);

            //Fill out non-static disorder Hamiltonians:
            nonstaticH[j] = xiNonstatic0*(h_r_x*p_x + h_r_y*p_y + h_r_z*p_z);
        }

        //Fill out ST matrices:
        for(int j = 0; j < circuitLength; j++){

            //Group single-site processes into pairs:
            MatrixX2cd tmpMat1 = 0.5*(staticH[j] + temporalH + nonstaticH[j]);
            MatrixX2cd tmpMat2 = 0.5*(staticH[(j+1)%circuitLength] + temporalH + nonstaticH[(j+1)%circuitLength]);
            MatrixX4cd tmpMat3 = Kronecker(tmpMat1,p_0) + Kronecker(p_0,tmpMat2);

            //Fill out Suzuki-Trotter matrices:
            if(j%2 != 0){
                stMatrices[j] = matrixExp(couplingTerm + tmpMat3, dt);
            }else{
                stMatrices[j] = matrixExp(couplingTerm + tmpMat3, dt/2);
            }
        }

        //Apply the Suzuki-Trotter matrices to the system:
        for(int k = 0; k < circuitLength; k+=2){
            QC1.applySites(k, stMatrices[k]);
        }
        for(int k = 1; k < circuitLength; k+=2){
            QC1.applySites(k, stMatrices[k]);
        }
        for(int k = 0; k < circuitLength; k+=2){
            QC1.applySites(k, stMatrices[k]);
        }

        //Apply measurement layer:
        if(measStrength0 != 0){
            for(int k = 0; k < circuitLength; k++){
                if(random_real(generator) < measProbability){
                    QC1.measureState(k);
                }
            }
        }

        //entropyTiming is in units of dt, so we just need to do an integer modulo check
        //Print the entropy:
       	if(i % entropyTiming == 0 && i*dt > startReadings){
      		QC1.findMI();
        }

        //Accumulating small errors due to computer precision requires regular renormalisation:
        if(i%100 == 0){
            QC1.normState();
        }
    }
    return 0;
}
