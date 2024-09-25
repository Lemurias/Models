#include <iostream>
#include <fstream>
#include <random>

#include "misc.hpp"

// The mod function in MATLAB returns the remainder after division of x by y.
uint32_t mod(uint32_t x, uint32_t y) {
	return ((x % y + y) % y);
}

// Function to generate normally distributed random numbers in an vector
std::vector<double> randn_vec(int size, double mean, double stddev)
{
    // Create a random number generator and normal distribution
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::normal_distribution<> d(mean, stddev); // Normal distribution with specified mean and stddev

    std::vector<double> numbers(size);

    for (int i = 0; i < size; ++i) {
        numbers[i] = d(gen); // Generate a random number from normal distribution
    }

    return numbers;
}

// Function to save double matrix NxM in a FILE
void saveMatrixInFile(const std::vector<std::vector<double>>& matrix, const std::string& fileName) {
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    uint32_t rows = static_cast<uint32_t>(matrix.size());
    uint32_t cols = static_cast<uint32_t>(matrix[0].size());

    /*    file << "Bitrate: " << param_Main.BitRate << "\n";
        file << "N samples per symbol: " << param_Main.n_samples_per_sym << "\n";
        file << "N samples per bit: " << param_Main.n_samples_per_sym << "\n";
        file << "N symbols: " << param_Main.n_symbols << "\n\n";
        */
    for (uint32_t i = 0; i < rows; ++i)
    {
        for (uint32_t j = 0; j < cols; ++j)
        {
            file << matrix[i][j];
            if (j < cols - 1) {
                file << ";";        // Separator between values ​​in the same row
            }
        }
        file << "\n";        // New line after each row
    }

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}

// Function to save double matrix NxM in a FILE
void saveMatrixInFile(const std::vector<double>& vect, const std::string& fileName) {
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    uint32_t cols = static_cast<uint32_t>(vect.size());

   for (uint32_t j = 0; j < cols; ++j)
   {
        file << vect[j] << "\n";
   }
   file << "\n";        // New line after each row

   file.close();
   std::cout << "Vector saved in " << fileName << std::endl;
}

// Function to save complex matrix NxM in a FILE
void saveMatrixInFile(const std::vector<std::vector<Complex>>& matrix, const std::string& fileName) {
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    uint32_t rows = static_cast<uint32_t>(matrix.size());
    uint32_t cols = static_cast<uint32_t>(matrix[0].size());

    /*    file << "Bitrate: " << param_Main.BitRate << "\n";
        file << "N samples per symbol: " << param_Main.n_samples_per_sym << "\n";
        file << "N samples per bit: " << param_Main.n_samples_per_sym << "\n";
        file << "N symbols: " << param_Main.n_symbols << "\n\n";
        */
    for (uint32_t i = 0; i < rows; ++i)
    {
        for (uint32_t j = 0; j < cols; ++j)
        {
            file << matrix[i][j];
            if (j < cols - 1) {
                file << ";";        // Separator between values ​​in the same row
            }
        }
        file << "\n";        // New line after each row
    }

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}

// Function to save complex matrix NxM in a FILE
void saveMatrixInFile(std::unique_ptr<std::unique_ptr<Complex[]>[]>& matrix, std::uint32_t rows, std::uint32_t cols, const std::string& fileName) {
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    /*    file << "Bitrate: " << param_Main.BitRate << "\n";
        file << "N samples per symbol: " << param_Main.n_samples_per_sym << "\n";
        file << "N samples per bit: " << param_Main.n_samples_per_sym << "\n";
        file << "N symbols: " << param_Main.n_symbols << "\n\n";
        */
    for (uint32_t i = 0; i < rows; ++i)
    {
        for (uint32_t j = 0; j < cols; ++j)
        {
            file << matrix[i][j];
            if (j < cols - 1) {
                file << ";";        // Separator between values ​​in the same row
            }
        }
        file << "\n";        // New line after each row
    }

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}


void saveMatrixInFile(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, std::uint32_t rows, std::uint32_t cols, const std::string& fileName)
{
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    /*    file << "Bitrate: " << param_Main.BitRate << "\n";
        file << "N samples per symbol: " << param_Main.n_samples_per_sym << "\n";
        file << "N samples per bit: " << param_Main.n_samples_per_sym << "\n";
        file << "N symbols: " << param_Main.n_symbols << "\n\n";
        */
    for (uint32_t i = 0; i < rows; ++i)
    {
        for (uint32_t j = 0; j < cols; ++j)
        {
            file << matrix[i][j];
            if (j < cols - 1) {
                file << ";";        // Separator between values ​​in the same row
            }
        }
        file << "\n";        // New line after each row
    }

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;

}

void saveMatrixInFile(std::unique_ptr<double[]>& vector, std::uint32_t cols, const std::string& fileName)
{
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    /*    file << "Bitrate: " << param_Main.BitRate << "\n";
        file << "N samples per symbol: " << param_Main.n_samples_per_sym << "\n";
        file << "N samples per bit: " << param_Main.n_samples_per_sym << "\n";
        file << "N symbols: " << param_Main.n_symbols << "\n\n";
        */
        for (uint32_t j = 0; j < cols; ++j)
        {
            file << vector[j];
            if (j < cols - 1) {
                file << ";";        // Separator between values ​​in the same row
            }
        }
        file << "\n";        // New line after each row
    
    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}

/***********************************************************************/


// Function to generate normally distributed random numbers in an m x n matrix
std::vector<std::vector<double>> randn(int m, int n, double mean, double stddev)
{
    // Initialize the random number generator with a seed based on the current time
    std::default_random_engine generator(static_cast<long unsigned int>(std::time(0)));
    std::normal_distribution<double> distribution(mean, stddev);

    // Create an m x n matrix filled with normally distributed random numbers
    std::vector<std::vector<double>> matrix(m, std::vector<double>(n));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = distribution(generator);
        }
    }

    return matrix;
}



// Function to compute the cumulative sum of a 2D array along rows
std::vector<std::vector<double>> cumsumRows(const std::vector<std::vector<double>>& A) {
    std::vector<std::vector<double>> result = A;
    for (size_t i = 0; i < A.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j) {
            sum += A[i][j];
            result[i][j] = sum;
        }
    }
    return result;
}

// Function to compute the cumulative sum of a 2D array along columns
std::vector<std::vector<double>> cumsumCols(const std::vector<std::vector<double>>& A) {
    if (A.empty()) return A;

    std::vector<std::vector<double>> result = A;
    for (size_t j = 0; j < A[0].size(); ++j) {
        double sum = 0.0;
        for (size_t i = 0; i < A.size(); ++i) {
            sum += A[i][j];
            result[i][j] = sum;
        }
    }
    return result;
}

std::vector<double> cumsum(const std::vector<double>& X) {
    std::vector<double> result(X.size());
    result[0] = X[0];

    for (size_t i = 1; i < X.size(); ++i) {
        result[i] = result[i - 1] + X[i];
    }

    return result;
}

// Function to perform element-wise multiplication of a scalar with a 2D array
std::vector<std::vector<double>> elementwise_multiply(const std::vector<std::vector<double>>& A, double scalar) {
    std::vector<std::vector<double>> result = A;  // Copy the original matrix

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            result[i][j] *= scalar;
        }
    }

    return result;
}

// Function to create an m x n matrix filled with ones
std::vector<std::vector<double>> ones(int m, int n) {
    // Create an m x n matrix initialized with 1.0
    return std::vector<std::vector<double>>(m, std::vector<double>(n, 1.0));
}

// Function to create an m x n matrix filled with zeros
std::vector<std::vector<double>> zeros(int m, int n) {
    // Create an m x n matrix initialized with 0.0
    return std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));
}

// Function to perform element-wise multiplication of two 2D vectors
std::vector<std::vector<double>> elementwise_multiply(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
        throw std::invalid_argument("Matrices must have the same dimensions.");
    }

    std::vector<std::vector<double>> result(mat1.size(), std::vector<double>(mat1[0].size()));
    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat1[i].size(); ++j) {
            result[i][j] = mat1[i][j] * mat2[i][j];
        }
    }
    return result;
}

// Function to compute the element-wise complex exponential of a vector
std::vector<Complex> complex_exponential(const std::vector<double>& vect) {
    std::vector<Complex> result(vect.size());
    for (size_t i = 0; i < vect.size(); ++i) {
        result[i] = std::exp(Complex(0, -2 * pi * vect[i])); // exp(-i * 2 * pi * mat[i][j])
    }
    return result;
}

// Function to compute the element-wise complex exponential of a 2D vector
std::vector<std::vector<Complex>> complex_exponential(const std::vector<std::vector<double>>& mat) {
    std::vector<std::vector<Complex>> result(mat.size(), std::vector<Complex>(mat[0].size()));
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[i].size(); ++j) {
            result[i][j] = std::exp(Complex(0, -2 * pi * mat[i][j])); // exp(-i * 2 * pi * mat[i][j])
        }
    }
    return result;
}

// Function to perform the final operation
std::vector<std::vector<Complex>> compute_E_laser(const std::vector<std::vector<Complex>>& E_l, const std::vector<std::vector<double>>& phi_tx, const std::vector<std::vector<double>>& t) {
    if (phi_tx.size() != t.size() || phi_tx[0].size() != t[0].size() || E_l.size() != t.size() || E_l[0].size() != t[0].size()) {
        throw std::invalid_argument("Matrices must have the same dimensions.");
    }

    std::vector<std::vector<double>> phi_tx_times_t = elementwise_multiply(phi_tx, t);
    std::vector<std::vector<Complex>> exp_term = complex_exponential(phi_tx_times_t);

    std::vector<std::vector<Complex>> E_laser(E_l.size(), std::vector<Complex>(E_l[0].size()));
    for (size_t i = 0; i < E_l.size(); ++i) {
        for (size_t j = 0; j < E_l[i].size(); ++j) {
            E_laser[i][j] = E_l[i][j] * exp_term[i][j];
        }
    }

    return E_laser;
}

// Function to create a 2D rotation matrix
std::vector<std::vector<double>> create_rotation_matrix(double theta_pol) {
    // Create a 2x2 matrix
    std::vector<std::vector<double>> rotation_matrix(2, std::vector<double>(2));

    // Fill the matrix with the appropriate trigonometric values
    rotation_matrix[0][0] = std::cos(theta_pol);
    rotation_matrix[0][1] = -std::sin(theta_pol);
    rotation_matrix[1][0] = std::sin(theta_pol);
    rotation_matrix[1][1] = std::cos(theta_pol);

    return rotation_matrix;
}


/***********************************************************************
    Functions to implement the FFT function
************************************************************************/
//
//unsigned int bitReverse(unsigned int x, int log2n) {
//    int n = 0;
//    int mask = 0x1;
//    for (int i = 0; i < log2n; i++) {
//        n <<= 1;
//        n |= (x & 1);
//        x >>= 1;
//    }
//    return n;
//}
//
//template<class Iter_T>
//void fft(Iter_T a, Iter_T b, int log2n)
//{
//    typedef typename iterator_traits<Iter_T>::value_type complex;
//    const complex J(0, 1);
//    int n = 1 << log2n;
//    for (unsigned int i = 0; i < n; ++i) {
//        b[bitReverse(i, log2n)] = a[i];
//    }
//    for (int s = 1; s <= log2n; ++s) {
//        int m = 1 << s;
//        int m2 = m >> 1;
//        complex w(1, 0);
//        complex wm = exp(-J * (pi / m2));
//        for (int j = 0; j < m2; ++j) {
//            for (int k = j; k < n; k += m) {
//                complex t = w * b[k + m2];
//                complex u = b[k];
//                b[k] = u + t;
//                b[k + m2] = u - t;
//            }
//            w *= wm;
//        }
//    }
//}
//
//
//// Another one optimized
//
//
//
//// Bit-reversal permutation
//void bit_reversal(CArray& a) {
//    size_t n = a.size();
//    size_t j = 0;
//    for (size_t i = 0; i < n; ++i) {
//        if (i < j) {
//            swap(a[i], a[j]);
//        }
//        size_t m = n;
//        while (j >= (m >>= 1)) {
//            j -= m;
//        }
//        j += m;
//    }
//}
//
//// Iterative FFT algorithm
//void fft(CArray& a) {
//    size_t n = a.size();
// //   assert((n & (n - 1)) == 0); // Ensure n is a power of 2
//
//    bit_reversal(a);
//
//    for (size_t len = 2; len <= n; len *= 2) {
//        double angle = -2 * pi / len;
//        Complex wlen(cos(angle), sin(angle));
//
//        for (size_t i = 0; i < n; i += len) {
//            Complex w(1);
//            for (size_t j = 0; j < len / 2; ++j) {
//                Complex u = a[i + j];
//                Complex t = w * a[i + j + len / 2];
//                a[i + j] = u + t;
//                a[i + j + len / 2] = u - t;
//                w *= wlen;
//            }
//        }
//    }
//}
//
//// Inverse FFT algorithm
//void ifft(CArray& a) {
//    size_t n = a.size();
//    // Conjugate the complex numbers
//    for (auto& x : a) {
//        x = conj(x);
//    }
//
//    // Apply FFT
//    fft(a);
//
//    // Conjugate the complex numbers again and normalize
//    for (auto& x : a) {
//        x = conj(x) / static_cast<double>(n);
//    }
//}


/*****************************************************/

std::vector<double> sign(std::vector<double>& matx) {


    std::vector<double> data(matx.size());

    for (size_t i = 0; i < matx.size(); ++i) {
        data[i] = (matx[i] >= 0) ? 1 : -1; // Apply the sign function
    }

    return data;


}


/*************************************************/
/************************************************************************/
// kron function in matlab
// Function to compute the Kronecker product
std::vector<std::vector<double>> kron(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    uint64_t a_rows = A.size();
    uint64_t a_cols = A[0].size();
    uint64_t b_rows = B.size();
    uint64_t b_cols = B[0].size();

    // Resulting matrix dimensions
    uint64_t r_rows = a_rows * b_rows;
    uint64_t r_cols = a_cols * b_cols;

    // Initialize the resulting matrix with zeros
    std::vector<std::vector<double>> result(r_rows, std::vector<double>(r_cols, 0));

    // Compute the Kronecker product
    for (uint64_t i = 0; i < a_rows; ++i) {
        for (uint64_t j = 0; j < a_cols; ++j) {
            for (uint64_t m = 0; m < b_rows; ++m) {
                for (uint64_t n = 0; n < b_cols; ++n) {
                    result[i * b_rows + m][j * b_cols + n] = A[i][j] * B[m][n];
                }
            }
        }
    }

    return result;
}



/****************************************************/
// Función recursiva FFT
void fft(ComplexVector& a) {
    size_t n = a.size();
    if (n <= 1) return;

    // Dividir en partes pares e impares
    ComplexVector a_even(n / 2);
    ComplexVector a_odd(n / 2);

    for (size_t i = 0; i < n / 2; ++i) {
        a_even[i] = a[i * 2];
        a_odd[i] = a[i * 2 + 1];
    }

    // Llamada recursiva
    fft(a_even);
    fft(a_odd);

    // Combinar
    for (size_t k = 0; k < n / 2; ++k) {
        Complex t = std::polar(1.0, -2 * pi * k / n) * a_odd[k];
        Complex u = a_even[k];
        a[k] = u + t;
        a[k + n / 2] = u - t;
    }
}

// Función recursiva IFFT
void ifft(ComplexVector& a) {
    size_t n = a.size();
    if (n <= 1) return;

    // Invertir el signo de los ángulos en la llamada recursiva
    ComplexVector a_even(n / 2);
    ComplexVector a_odd(n / 2);

    for (size_t i = 0; i < n / 2; ++i) {
        a_even[i] = a[i * 2];
        a_odd[i] = a[i * 2 + 1];
    }

    // Llamada recursiva
    ifft(a_even);
    ifft(a_odd);

    // Combinar, pero ahora el ángulo es positivo
    for (size_t k = 0; k < n / 2; ++k) {
        Complex t = std::polar(1.0, 2 * pi * k / n) * a_odd[k];
        Complex u = a_even[k];
        a[k] = u + t;
        a[k + n / 2] = u - t;
    }

    // Normalizar los resultados
    for (size_t i = 0; i < n; ++i) {
        a[i] /= static_cast<double>(n);
    }
}


// Función de FFT iterativa (Cooley-Tukey)
void fft_o(ComplexVector& a) {
    size_t n = a.size();
    if (n <= 1) return;

    // Reordenar los datos en orden bit-reversal
    size_t log_n = static_cast<size_t>(std::log2(n));
    for (size_t i = 0; i < n; ++i) {
        size_t j = 0;
        for (size_t k = 0; k < log_n; ++k) {
            j |= ((i >> k) & 1) << (log_n - 1 - k);
        }
        if (j > i) std::swap(a[i], a[j]);
    }

    // Aplicar el algoritmo FFT
    for (size_t len = 2; len <= n; len *= 2) {
        double angle = -2 * pi / len;
        Complex wlen(std::cos(angle), std::sin(angle));
        for (size_t i = 0; i < n; i += len) {
            Complex w(1);
            for (size_t j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex t = w * a[i + j + len / 2];
                a[i + j] = u + t;
                a[i + j + len / 2] = u - t;
                w *= wlen;
            }
        }
    }
}

// Función de IFFT iterativa (reciprocal de FFT)
void ifft_o(ComplexVector& a) {
    size_t n = a.size();
    if (n <= 1) return;

    // Tomar la FFT
    fft(a);

    // Normalizar y ajustar la fase
    for (size_t i = 0; i < n; ++i) {
        a[i] /= static_cast<double>(n);
    }

    // Invertir la fase de los elementos
    for (size_t i = 0; i < n; ++i) {
        a[i] = std::conj(a[i]);
    }
}

// Función para imprimir un vector de std::complex<double>
void printArray(const ComplexVector& arr) {
    for (const auto& val : arr) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}


// Function to convert a vector of double to a vector of std::complex<double>
ComplexVector convertToComplex(const DoubleVector& realVector) {
    ComplexVector complexVector;
    complexVector.reserve(realVector.size()); // Reserve space to avoid relocation

    for (double value : realVector) {
        complexVector.emplace_back(value, 0.0); // Convert to std::complex<double>
    }

    return complexVector;
}

// Function to print a vector of std::complex<double>
void printComplexVector(const ComplexVector& vec) {
    for (const auto& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


// multiply two double matrix
DoubleMatrix multiply_d_matrices(const DoubleMatrix& A, const DoubleMatrix& B) {
    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t rows_B = B.size();
    size_t cols_B = B[0].size();

    // check if it´s posible to multiply 
    if (cols_A != rows_B) {
        throw std::invalid_argument("Number of columns in A must be equal to the number of rows in B.");
    }

    // Create the result matrix with appropriate dimensions
    DoubleMatrix C(rows_A, std::vector<double>(cols_B, (0.0)));

    // Perform matrix multiplication
    for (size_t i = 0; i < rows_A; ++i) {
        for (size_t j = 0; j < cols_B; ++j) {
            for (size_t k = 0; k < cols_A; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// multiply two complex matrix
ComplexMatrix multiply_matrices(const ComplexMatrix& A, const ComplexMatrix& B) {
    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t rows_B = B.size();
    size_t cols_B = B[0].size();

    // check if it´s posible to multiply 
    if (cols_A != rows_B) {
        throw std::invalid_argument("Number of columns in A must be equal to the number of rows in B.");
    }

    // Create the result matrix with appropriate dimensions
    ComplexMatrix C(rows_A, std::vector<Complex>(cols_B, Complex(0.0, 0.0)));

    // Perform matrix multiplication
    for (size_t i = 0; i < rows_A; ++i) {
        for (size_t j = 0; j < cols_B; ++j) {
            for (size_t k = 0; k < cols_A; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// multiply a double matrix by a complex matrix
ComplexMatrix multiply_matrices(const DoubleMatrix& A, const ComplexMatrix& B) {
    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t rows_B = B.size();
    size_t cols_B = B[0].size();

    // check if it´s posible to multiply 
    if (cols_A != rows_B) {
        throw std::invalid_argument("Number of columns in A must be equal to the number of rows in B.");
    }

    // Create the result matrix with appropriate dimensions
    ComplexMatrix C(rows_A, std::vector<Complex>(cols_B, Complex(0.0, 0.0)));

    // Perform matrix multiplication
    for (size_t i = 0; i < rows_A; ++i) {
        for (size_t j = 0; j < cols_B; ++j) {
            for (size_t k = 0; k < cols_A; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}




// built a 2xN matrix from a double complex vector
ComplexMatrix create_matrix(const ComplexVector& vec) {
    size_t N = vec.size();
    ComplexMatrix matrix(2, std::vector<Complex>(N, Complex(0.0, 0.0)));

    // Fill the first row with the values ​​of the vector
    for (size_t i = 0; i < N; ++i) {
        matrix[0][i] = vec[i];
    }

    // The second row is already filled with zeros by default

    return matrix;
}

// convert a double matrix to complex double matrix
ComplexMatrix convert_to_complex_matrix(const DoubleMatrix& double_matrix) {
    size_t N = double_matrix.size();
    size_t M = double_matrix[0].size();
    ComplexMatrix complex_matrix(N, std::vector<Complex>(M, Complex(0.0, 0.0)));

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            complex_matrix[i][j] = Complex(double_matrix[i][j], 0.0);
        }
    }

    return complex_matrix;
}

/********************************************************/
// Function to compute the FFTShift function

std::vector<double> fftshift(const std::vector<double>& input) {
    std::vector<double> output(input);
    uint32_t N = output.size();

    // Calcula el punto medio
    uint32_t mid = (N + 1) / 2;

    // Realiza la rotación circular para mover el cero a la mitad del vector
    std::rotate(output.begin(), output.begin() + mid, output.end());

    return output;
}


void fft_shift(std::vector<double>& input) {
    uint32_t n = input.size();
    uint32_t half_n = n / 2;

    if (n % 2 == 0) {
        std::rotate(input.begin(), input.begin() + half_n, input.end());
    }
    else {
        std::rotate(input.begin(), input.begin() + half_n + 1, input.end());
    }
    return;
}

void fft_shift(std::vector<Complex>& input) {
    uint32_t n = input.size();
    uint32_t half_n = n / 2;

    if (n % 2 == 0) {
        std::rotate(input.begin(), input.begin() + half_n, input.end());
    }
    else {
        std::rotate(input.begin(), input.begin() + half_n + 1, input.end());
    }
    return;
}

/********************************************************/
// Function to compute the FFTShift function
//void fftshift(double *data, uint32_t n) {
//    // event number of elements
//    if (n % 2 == 0)
//        std::rotate(&data[0], &data[n >> 1], &data[n]);
//    // odd number of elements
//    else
//        std::rotate(&data[0], &data[(n >> 1)+1], &data[n]);
//
//}

// Optimized function to perform ifftshift on a 1D array
std::vector<double> ifftshift(const std::vector<double>& data) {
    int n = data.size();
    int shift_index = (n + 1) / 2;  // Shift index for ifftshift
    
    std::vector<double> shifted_data;
    shifted_data.reserve(n); // Reserve memory to avoid reallocations
    
    // Append second half of the data
    shifted_data.insert(shifted_data.end(), data.begin() + shift_index, data.end());
    // Append first half of the data
    shifted_data.insert(shifted_data.end(), data.begin(), data.begin() + shift_index);
    
    return shifted_data;
}

// Optimized function to perform ifftshift on a 1D array
std::vector<Complex> ifftshift(const std::vector<Complex>& data) {
    int n = data.size();
    int shift_index = (n + 1) / 2;  // Shift index for ifftshift
    
    std::vector<Complex> shifted_data;
    shifted_data.reserve(n); // Reserve memory to avoid reallocations
    
    // Append second half of the data
    shifted_data.insert(shifted_data.end(), data.begin() + shift_index, data.end());
    // Append first half of the data
    shifted_data.insert(shifted_data.end(), data.begin(), data.begin() + shift_index);
    
    return shifted_data;
}

// Optimized function to perform ifftshift on a 2D matrix
std::vector<std::vector<double>> ifftshift(const std::vector<std::vector<double>>& data) {
    int rows = data.size();
    int cols = data[0].size();
    
    int row_shift = (rows + 1) / 2;  // Shift index for rows
    int col_shift = (cols + 1) / 2;  // Shift index for columns
    
    std::vector<std::vector<double>> shifted_data(rows, std::vector<double>(cols));
    
    // Rebuild the matrix in an optimized way
    for (int i = 0; i < rows; ++i) {
        // Calculate the new row index after shifting
        int new_row = (i + row_shift) % rows;
        
        // Copy the second half of the columns
        std::copy(data[i].begin() + col_shift, data[i].end(), shifted_data[new_row].begin());
        // Copy the first half of the columns
        std::copy(data[i].begin(), data[i].begin() + col_shift, shifted_data[new_row].begin() + (cols - col_shift));
    }

    return shifted_data;
}

// Function to perform element-wise division of a 2D matrix by a scalar
std::vector<std::vector<double>> elementWiseDivision(const std::vector<std::vector<double>>& matrix, double scalar) {
    // Get the number of rows and columns in the matrix
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // Create a result matrix with the same dimensions
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

    // Perform the element-wise division
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i][j] = matrix[i][j] / scalar;
        }
    }

    return result;
}

double variance(const std::vector<double>& data, bool isSample){
    int n = data.size();
    if (n == 0) return 0.0;  // Avoid division by zero
    
    // Step 1: Calculate the mean
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;

    // Step 2: Compute the variance
    double var = 0.0;
    for (const double& x : data) {
        var += std::pow(x - mean, 2);
    }

    // If isSample is true, use n-1 for sample variance, else use n for population variance
    return var / (isSample ? n - 1 : n);
}