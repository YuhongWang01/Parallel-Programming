#include <vector>
#include <iostream>
#include<semaphore.h>
#include <emmintrin.h>
#include <immintrin.h>
using namespace std;

class Grobner_Matrix {
    /*
     * Grobner_Matrixʵ���˶�Grobner����������жԾ����ѹ��
     * ԭ����Ϊn*n����Ϊֻ�ڣ�0,1�����ϼ��㣬���Կ���ѹ����n*n/32��С�ľ���.
     */
public:
    int n;//������
    int m;//������

    int m_;//ѹ���������
    vector<int> row_index;//������,���ڼ�¼��Ԫ�ӻ���Ԫ���Ƿ����
    int** matrix;//ԭ����
    Grobner_Matrix();
    explicit Grobner_Matrix(int n, int m);
    ~Grobner_Matrix();
    void init();
    void set_bit(int i, int j);
    int xor_line(Grobner_Matrix&, int i, int j);
    int Simd_xor_line(Grobner_Matrix&, int i, int j);
    int get_max_bit(int i);
    void input_line(int i, vector<int>&);
    int dxor(vector<int>l, int j);
    vector<int> get_5_line(int s);
    vector<int> get_line(int s, int num);
    void print_line(int i);

    void clear();
};

Grobner_Matrix::Grobner_Matrix() {
    n = 0;
    m = 0;
}

Grobner_Matrix::Grobner_Matrix(int n, int m) {
    this->n = n;
    this->m = m;
    init();
}

Grobner_Matrix::~Grobner_Matrix() {
    delete[] matrix;
}

void Grobner_Matrix::init() {
    matrix = new int* [n];
    m_ = (m - 1) / 32 + 1;

    for (int i = 0; i < n; i++) {
        matrix[i] = new int[m_];
    }
    row_index.resize(n, -1);
}

void Grobner_Matrix::set_bit(int i, int j) {
    //����i�еĵ�j��bit��1
    matrix[i][j / 32] |= (1 << (j % 32));
}

int Grobner_Matrix::xor_line(Grobner_Matrix& grobnerMatrix, int i, int j) {
    //�����������ķ���λ
    int max_bit = -1;
    for (int k = 0; k < m_; k++) {
        matrix[j][k] ^= grobnerMatrix.matrix[i][k];
    }
    max_bit = get_max_bit(j);
    return max_bit;
}
int Grobner_Matrix::dxor(vector<int> l, int j) {
    int max_bit = 0;
    for (int k = 0; k < m_; k++) {
        matrix[j][k] ^= l[k];
    }
    max_bit = get_max_bit(j);
    return max_bit;
}

void Grobner_Matrix::input_line(int i, vector<int>& line) {
    //����i�е��������뵽������
    for (int j : line) {
        set_bit(i, j);
    }
    row_index[i] = i;
}

vector<int> Grobner_Matrix::get_5_line(int s) {
    //��s��ʼ��ȡ��ǰ����
    vector<int> top5_line;
    top5_line.resize(row_index.size() - 5 * s >= 5 ? 5 : row_index.size() - 5 * s);
    for (int i = 0; i < top5_line.size(); i++) {
        top5_line[i] = row_index[i + 5 * s];
    }
    return top5_line;
}

vector<int> Grobner_Matrix::get_line(int s, int num) {
    //��s��ʼ��ȡ��ǰnum��
    vector<int> top_line;
    top_line.resize(row_index.size() - num * s >= num ? num : row_index.size() - num * s);
    for (int i = 0; i < top_line.size(); i++) {
        top_line[i] = row_index[i + num * s];
    }
    return top_line;
}

int Grobner_Matrix::get_max_bit(int i) {
    //��ȡ��i�����ķ���λ
    int max_bit = -1;
    for (int j = 0; j < m_; j++) {
        for (int k = 0; k < 32; k++) {
            if ((matrix[i][j] & (1 << k)) != 0)
                max_bit = j * 32 + k;
        }
    }
    return max_bit;
}

void Grobner_Matrix::print_line(int i) {
    //����i�и���Ϊ1��λ���Ӵ�С���
    for (int j = 0; j < m_; j++) {
        if (matrix[i][j] != 0) {
            for (int k = 0; k < 32; k++) {
                if (matrix[i][j] & (1 << k))
                    cout << j * 32 + k << " ";
            }
        }
    }
    cout << endl;
}

void Grobner_Matrix::clear() {
    //��վ���
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m_; j++) {
            matrix[i][j] = 0;
        }
    }
    row_index.clear();
}

int Grobner_Matrix::Simd_xor_line(Grobner_Matrix& grobnerMatrix, int i, int j) {
    //�����������ķ���λ
    int max_bit = 0;
    int k = 0;
    for (k = 0; k + 4 <= m_; k += 4) {
        __m128i mjk = _mm_load_si128((__m128i*) & matrix[j][k]);
        __m128i mik = _mm_load_si128((__m128i*) & grobnerMatrix.matrix[i][k]);
        mjk = _mm_xor_si128(mjk, mik);
        _mm_store_si128((__m128i*) & matrix[j][k], mjk);
        //matrix[j][k] ^= grobnerMatrix.matrix[i][k];
    }
    while (k < m_) {
        matrix[j][k] ^= grobnerMatrix.matrix[i][k];
        k++;
    }
    max_bit = get_max_bit(j);
    return max_bit;
}