#include "PicReader.h"
#include "huffmancode.h"
#include<iostream>
#include <fstream>
#include <stdio.h>
#include<cmath>
#include<map>
#include<bitset>
using namespace std;

class jpegcompress {
    private:
        /*BYTE*&, ����BYTEָ������,���ڽ���������Ϣ
        * img_width ͼ��Ŀ����Ϣ
        * img_height ͼ��ĸ߶���Ϣ
          * ������Ϣÿ�ĸ�һ��(R G B A)    */
        PicReader imread;
        BYTE* data = nullptr;
        UINT img_width, img_height;
        
        //������������ұ�
        map<unsigned char, string> mp_DC_lumin;  //����ֱ��
        map<unsigned char, string> mp_AC_lumin;  //���Ƚ���
        map<unsigned char, string> mp_DC_chro;   //ɫ��ֱ��
        map<unsigned char, string> mp_AC_chro;   //ɫ�Ƚ���

    public:
        string tmp = "";  //���д������ʱ��8λ01��

        //����ͼƬ
        void read_picture(char* argv[]);
        //�ļ�ͷ��д��
        void prewrite(ofstream &fout);
        //��������huffman������ұ�
        void huffmancode();
        //����huffman������ұ�
        void set_huffmancode_list(int chooselist, const char* code_cnt, const unsigned char* code_value);
        //д��ѹ�����ͼ����Ϣ
        void writecode(int zz[64], int chooselist, ofstream &fout);
        //�õ���8*8С���01������
        void getcode(int zz[64], int chooselist, string& result);
        //����ԭͼ������Ϣ
        void process_pixel(ofstream &fout);
        //д���ļ��������ͷſռ�
        void end(ofstream &fout);
};

int zz_Y[N * N] = { 0 }, zz_Cb[N * N] = { 0 }, zz_Cr[N * N] = { 0 };  //zigzagɨ������������
int preY = 0, preCb = 0, preCr = 0;  //�����һС���һ�����ݵ�ԭֵ�����ڲ�֣�

//����ͼƬ
void jpegcompress::read_picture(char* argv[])
{
    imread.readPic(argv[2]);
    imread.getData(data, img_width, img_height);
}

//����ÿС��RLE����ĵڶ�������Ӧ�ı���
string getcode_second(int val)
{
    string codeval = "";
    for (int k = abs(val); k > 0; k = k / 2) {
        codeval = codeval + (k % 2 == 1 ? '1' : '0');
    }

    if (val < 0)  //�����������
        for (int k = 0; k < (int)codeval.length(); k++)
            if (codeval[k] == '1')
                codeval[k] = '0';
            else
                codeval[k] = '1';

    //�ַ�������
    for (int i = 0; i < (int)codeval.length() / 2; i++) {
        char ch = codeval[i];
        codeval[i] = codeval[codeval.length() - i - 1];
        codeval[codeval.length() - i - 1] = ch;
    }
    return codeval;
}

static string to_binary_str(int code, int bitslen, string buf) {
    int mask = 1 << (bitslen - 1);
    for (int i = 0; i < bitslen; i++) {
        if (code & mask)
            buf += '1';
        else
            buf += '0';
        mask >>= 1;
    }
    buf[bitslen] = '0';
    return buf;
}

//����huffman������ұ�
//chooselist����ֱ����������ɫ�ȣ�code_cntָ��ͬ���ȱ���ĸ������飬code_valueΪ��������ַ�
void jpegcompress::set_huffmancode_list(int chooselist, const char* code_cnt, const unsigned char* code_value)
{
    int sym = 0, code = 0;
    string buf = "";
    for (int i = 0; i < 16; i++) {
        int num = (int)code_cnt[i];
        int bitslen = i + 1;
        for (int j = 0; j < num; j++) {
            unsigned char symbol = code_value[sym];
            string result = to_binary_str(code, bitslen, buf);

            if (chooselist == 1) 
                mp_DC_lumin.insert(make_pair(symbol, result));
            else if (chooselist == 2)
                mp_AC_lumin.insert(make_pair(symbol, result));
            else if (chooselist == 3)  
                mp_DC_chro.insert(make_pair(symbol, result));
            else  
                mp_AC_chro.insert(make_pair(symbol, result));
      
            code++;
            sym++;
        }
        code <<= 1;
    }
}

//��������huffman������ұ�
void jpegcompress::huffmancode()
{
    set_huffmancode_list(1, DC_lumin, DC_lumin_value);
    set_huffmancode_list(2, AC_lumin, AC_lumin_value);
    set_huffmancode_list(3, DC_chro, DC_chro_value);
    set_huffmancode_list(4, AC_chro, AC_chro_value);
}

//�õ���8*8С���01������
//zzΪ��zigzag��ʽ���гɵ�һά���飬chooselist����ֱ����������ɫ�ȣ�result�������
void jpegcompress::getcode(int zz[64], int chooselist, string& result)
{
    result = "";
    //numof0��ʾ��С�����һ����ǰ��0�ĸ�����val��ʾ�ö�ĩβ����
    int numof0 = 0, val = 0;  
    int end = 64;  //end��ʾEOB��λ��
    for (int i = 63; i >= 0; i--)
        if (zz[i] != 0) {
            end = i + 1;
            break;
        }

    //DC����
    string codeval = getcode_second(zz[0]);
    int len = codeval.length();
    if (len != 0) {
        if (chooselist == 1)
            result += mp_DC_lumin[(unsigned char)len];
        else
            result += mp_DC_chro[(unsigned char)len];
    }
    else {
        if (chooselist == 1)
            result += mp_DC_lumin[(unsigned char)0x00];
        else
            result += mp_DC_chro[(unsigned char)0x00];
    }
    result += codeval;

    //AC����
    for (int i = 1; i < end; i++) {
        if (zz[i] == 0 && numof0 < 15)  //ÿС���Է�0����β�������16����
            numof0++;
        else {
            val = zz[i];
            string codeval = getcode_second(val);  //ÿС��RLE����ĵڶ�������Ӧ�ı���
            int len = codeval.length();  //BIT������м�һ����
            

            if (chooselist == 1) {
                    result += mp_AC_lumin[(unsigned char)((numof0 << 4) + len)];
            }
            else {
                    result += mp_AC_chro[(unsigned char)((numof0 << 4) + len)];
            }
            result += codeval;

            numof0 = 0;
        }
    }
}

//д������������
void print_huffman(char tag, const char  valuelen[], const unsigned char value[], ofstream &fout)
{
    fout.put(tag);  //��ID�ͱ�����
    int sum = 0;    //����value��Ԫ�ظ���Ϊvaluelen����Ԫ��֮��
    for (int i = 0; i < 16; i++) {
        fout.put(valuelen[i]);
        sum += (int)valuelen[i];
    }
    for (int i = 0; i < sum; i++)
        fout.put(value[i]);
}

//���������������
void print_Q(const unsigned char Q[N][N], ofstream &fout)
{
    int tx = 0, ty = 0, dx = -1, dy = 1, cnt = 0;
    while (tx != N - 1 || ty != N - 1) {
        fout.put(Q[tx][ty]);
        if ((tx == 0 || tx == N - 1) && ty % 2 == 0) {
            ty++;
            dx *= -1;
            dy *= -1;
            continue;
        }
        if ((ty == 0 || ty == N - 1) && tx % 2 == 1) {
            tx++;
            dx *= -1;
            dy *= -1;
            continue;
        }
        tx += dx;
        ty += dy;
    }
    fout.put(Q[N - 1][N - 1]);
}

//�ļ�ͷ��д��
void jpegcompress::prewrite(ofstream &fout)
{
    /*SOI ͼ��ʼ*/
    fout.put((char)0xff), fout.put((char)0xd8);
    fout.put((char)0xff), fout.put((char)0xe0);  //APP0
    fout.put(0x00), fout.put(0x10);
    fout.put(0x4a), fout.put(0x46), fout.put(0x49), fout.put(0x46), fout.put(0x00);
    fout.put(0x01), fout.put(0x01), fout.put(0x00), fout.put(0x00), fout.put(0x01), fout.put(0x00), fout.put(0x01);
    fout.put(0x00), fout.put(0x00);

    /*DQT ������*/
    fout.put((char)0xff), fout.put((char)0xdb);  //��Ǵ���
    fout.put((char)0x00), fout.put((char)0x84);  //���ݳ���(�������ֶ�)
    fout.put(0x00);  //���ȼ�������ID
    print_Q(Qy, fout);    //���������������Qy
    fout.put(0x01);  //���ȼ�������ID
    print_Q(Qc, fout);    //���������������Qc

    /*SOF0 ֡ͼ��ʼ*/
    fout.put((char)0xff), fout.put((char)0xc0);  //��Ǵ���
    fout.put(0x00), fout.put(0x11), fout.put(0x08);
    fout.put(0x02); fout.put(0x00); //fout << img_height;  //ͼ��߶�(2�ֽ�)
    fout.put(0x02); fout.put(0x00); //fout << img_width;  //ͼ����(2�ֽ�)
    fout.put(0x03);  //��ɫ������
    //��ɫ������Ϣ
    fout.put(0x01), fout.put(0x11), fout.put(0x00); 
    fout.put(0x02), fout.put(0x11), fout.put(0x01);
    fout.put(0x03), fout.put(0x11), fout.put(0x01);

    /*DHT �����������*/
    fout.put((char)0xff), fout.put((char)0xc4);  //��Ǵ���
    fout.put((char)0x01), fout.put((char)0xa2);  //���ݳ���
    print_huffman(0x00, DC_lumin, DC_lumin_value, fout);
    print_huffman(0x10, AC_lumin, AC_lumin_value, fout);
    print_huffman(0x01, DC_chro, DC_chro_value, fout);
    print_huffman(0x11, AC_chro, AC_chro_value, fout);

    /*SOS ɨ�迪ʼ*/
    fout.put((char)0xff), fout.put((char)0xda);  //��Ǵ���
    fout.put(0x00), fout.put(0x0c);  //���ݳ���
    fout.put(0x03);  //��ɫ������
    fout.put(0x01), fout.put(0x00), fout.put(0x02), fout.put(0x11), fout.put(0x03), fout.put(0x11);  //��ɫ������Ϣ (?)
    fout.put(0x00), fout.put(0x3f), fout.put(0x00);  //�̶�ֵ
}

//д��ѹ�����ͼ����Ϣ
//zzΪ��zigzag��ʽ���гɵ�һά���飬chooselist����ֱ����������ɫ��
void jpegcompress::writecode(int zz[64], int chooselist, ofstream &fout)
{
    string result = "";
    getcode(zz, chooselist, result);  //��ȡEOBǰ�ı��룬�����ַ���result��
    if (chooselist == 1)
        result += mp_AC_lumin[(unsigned char)0x00];
    else
        result += mp_AC_chro[(unsigned char)0x00];
    //cout << result;
    //��result�����ֽ�д���ļ�
    for (int pos = 0; pos < (int)result.length(); pos++) {
        tmp += result[pos];
        if ((int)tmp.length() >= 8) {   //tmp��8λ�򽫸��ֽڱ���תΪ�ַ����
            //cout << n << endl;
            bitset<8> n(tmp);
            unsigned char ch = (unsigned char)(n.to_ulong());
            fout.put(ch);
            if (ch == 0xff)
                fout.put(0x00);
            tmp = "";   //���tmp
        }
    }
}

double alpha(int x)
{
    if (x == 0)
        return 1.0 / sqrt(8);
    else
        return 0.5;
}

//����ԭͼ������Ϣ
void jpegcompress::process_pixel(ofstream &fout)
{
    int count = 0;  //��ǰ��count��8*8С��
    for (DWORD i = 0; i < img_height; i += 8)
        for (DWORD j = 0; j < img_width; j += 8) {
            //ͼ��ָ��ÿ��8*8С��
            count++;
            double Y[N][N] = { 0 }, Cb[N][N] = { 0 }, Cr[N][N] = { 0 };
            for (DWORD u = 0; u < 8; u++)
                for (DWORD v = 0; v < 8; v++) {
                    int k = (i + u) * img_width * 4 + (j + v) * 4;
                    //��ɫ�ռ�ת��RGB->YCbCr
                    int R = data[k], G = data[k + 1], B = data[k + 2], A = data[k + 3];
                    Y[u][v] = 0.29871 * R + 0.58661 * G + 0.11448 * B - 128;
                    Cb[u][v] = -0.16874 * R - 0.33126 * G + 0.5 * B;
                    Cr[u][v] = 0.5 * R - 0.41869 * G - 0.08131 * B;
                }

            //��ɢ���ұ任
            double f_Y[N][N] = { 0 }, f_Cb[N][N] = { 0 }, f_Cr[N][N] = { 0 };
            for (int u = 0; u < N; u++)
                for (int v = 0; v < N; v++) {
                    double tmp1 = 0, tmp2 = 0, tmp3 = 0;
                    for (int tu = 0; tu < N; tu++)
                        for (int tv = 0; tv < N; tv++) {
                            tmp1 += Y[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                            tmp2 += Cb[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                            tmp3 += Cr[tu][tv] * cos(1.0*(2 * tu + 1) * u * Pi / 16.0) * cos(1.0*(2 * tv + 1) * v * Pi / 16.0);
                        }
                    f_Y[u][v] = alpha(u) * alpha(v) * tmp1;
                    f_Cb[u][v] = alpha(u) * alpha(v) * tmp2;
                    f_Cr[u][v] = alpha(u) * alpha(v) * tmp3;
                }

            //��������
            int qua_Y[N][N], qua_Cb[N][N], qua_Cr[N][N];
            for (int u = 0; u < N; u++)
                for (int v = 0; v < N; v++) {
                    qua_Y[u][v] = (int)round(f_Y[u][v] / (int)Qy[u][v]);
                    qua_Cb[u][v] = (int)round(f_Cb[u][v] / (int)Qc[u][v]);
                    qua_Cr[u][v] = (int)round(f_Cr[u][v] / (int)Qc[u][v]);
                }
            
            //Zigzagɨ������ (��������Ķ�ά����ת���һ��һά����)

            int tmpY = qua_Y[0][0];
            int tmpCb = qua_Cb[0][0];
            int tmpCr = qua_Cr[0][0];
            zz_Y[0] = qua_Y[0][0] - preY;
            zz_Cb[0] = qua_Cb[0][0] - preCb;
            zz_Cr[0] = qua_Cr[0][0] - preCr;
            preY = tmpY;
            preCb = tmpCb;
            preCr = tmpCr;

            int tx = 0, ty = 1, dx = 1, dy = -1, cnt = 1;
            while (tx != N - 1 || ty != N - 1) {
                zz_Y[cnt] = qua_Y[tx][ty];
                zz_Cb[cnt] = qua_Cb[tx][ty];
                zz_Cr[cnt] = qua_Cr[tx][ty];
                cnt++;
                if ((tx == 0 || tx == N - 1) && ty % 2 == 0) {
                    ty++;
                    dx *= -1;
                    dy *= -1;
                    continue;
                }
                if ((ty == 0 || ty == N - 1) && tx % 2 == 1) {
                    tx++;
                    dx *= -1;
                    dy *= -1;
                    continue;
                }
                tx += dx;
                ty += dy;
            }
            zz_Y[cnt] = qua_Y[N - 1][N - 1];
            zz_Cb[cnt] = qua_Cb[N - 1][N - 1];
            zz_Cr[cnt] = qua_Cr[N - 1][N - 1];

            //д��ѹ�����ͼ����Ϣ
            writecode(zz_Y, 1, fout);
            writecode(zz_Cb, 2, fout);
            writecode(zz_Cr, 2, fout);
        }
        if (tmp != "") {    //������ʣ�಻��8λ��
            bitset<8> n(tmp);
            fout.put((unsigned char)n.to_ulong());
        }
}

//д���ļ��������ͷſռ�
void jpegcompress::end(ofstream &fout)
{
    fout.put((char)0xff);
    fout.put((char)0xd9);
    delete[] data;
    data = nullptr;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Please make sure the number of parameters is correct." << endl;
        return -1;
    }
    if (strcmp(argv[1], "-read") && strcmp(argv[1], "-compress")) {
        cerr << "Unknown parameter!\nCommand list:\nzip/unzip" << endl;
        return -1;
    }

    jpegcompress jpeg;

    if (strcmp(argv[1], "-compress") == 0) {
        //�ļ���
        ifstream fin(argv[2], ios::binary); // �Զ����Ʒ�ʽ�������ļ�
        if (!fin) {  // ���������Ϣ���˳�
            cerr << "Can not open the input file!" << endl;
            return -1;
        }

        ofstream fout("lena.jpg", ios::binary);
        if (!fout) {
            cerr << "Can not open the output file!" << endl;
            return -1;
        }

        jpeg.read_picture(argv);
        jpeg.prewrite(fout);
        jpeg.huffmancode();
        jpeg.process_pixel(fout);
        jpeg.end(fout);
        cout << "complete!" << endl;

        fin.close();
        fout.close();
    }
    else if (strcmp(argv[1], "-read") == 0) {
        PicReader imread;
        BYTE* data = nullptr;
        UINT x, y;

        imread.readPic(argv[2]);
        imread.getData(data, x, y);
        imread.showPic(data, x, y);
        delete[] data;
    }
 
    return 0;
}